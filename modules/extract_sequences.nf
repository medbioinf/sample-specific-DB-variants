// Extracts coding sequences (CDS) from the reference genome and from a sample-specific genome containing protein-altering mutations.
// Requires 5 items:
//      reference genome: .fa + .fa.fai
//      BED file with CDS regions of affected transcripts (extract_cds_bed.nf output): .bed 
//      protein altering variants (protein_altering_bundle output): .vcf.gz + .vcf.gz.tbi

process extract_sequences {
    tag "$sample_id"
    publishDir "${params.outdir}/06_vcf/05_extract_sequences", mode: 'copy', overwrite: true, dynamic: true
    cpus { params.threads ?: 8 }
    container (params.toolbox_image ?: 'biocontainers/bcftools:v1.20-1--deb_cv1')
    
    input:
    tuple val(sample_id),
       path(reference_fasta),
       path(reference_fai),
       path(cds_affected_bed),
       path(filtered_vcf),
       path(filtered_vcf_index)

    when:
    reference_fasta.exists() && reference_fai.exists() &&
    cds_affected_bed.exists() &&
    filtered_vcf.exists() && filtered_vcf_index.exists()

    output:
    path("${sample_id}.original_cds.fa"), emit: original_cds_fa
    path("genome_mutated.fa"),            emit: genome_mutated_fa
    path("${sample_id}.mutated_cds.fa"),  emit: mutated_cds_fa

    script:
    """
    set -euo pipefail

    frame_bed="${sample_id}.cds_affected.frame.bed"
    python3 - <<'PY' "${cds_affected_bed}" "\${frame_bed}"
import sys
from collections import defaultdict
from pathlib import Path

source = Path(sys.argv[1])
out = Path(sys.argv[2])

by_tx = defaultdict(list)
with source.open() as handle:
    for line in handle:
        if not line.strip():
            continue
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) < 6:
            continue
        chrom, start, end, tx = parts[:4]
        strand = parts[5] if len(parts) > 5 and parts[5] else "+"
        try:
            phase = int(parts[6]) if len(parts) > 6 and parts[6] not in {".",""} else 0
        except Exception:
            phase = 0
        try:
            s = int(start)
            e = int(end)
        except ValueError:
            continue
        by_tx[tx].append((s, e, chrom, strand, phase))

with out.open("w") as handle:
    for tx, entries in by_tx.items():
        strand = entries[0][3]
        sorted_entries = sorted(entries, key=lambda x: x[0], reverse=(strand == "-"))
        for s, e, chrom, strand, phase in sorted_entries:
            # Apply phase to stay in-frame: plus strand -> start+phase; minus strand -> end-phase
            if strand == "-":
                adj_start = s
                adj_end = e - phase
            else:
                adj_start = s + phase
                adj_end = e
            if adj_end <= adj_start:
                continue
            handle.write(f"{chrom}\\t{adj_start}\\t{adj_end}\\t{tx}\\t0\\t{strand}\\n")
PY

    # Step 1: extract original CDS from reference genome
    bedtools getfasta \\
        -fi ${reference_fasta} \\
        -bed \${frame_bed} \\
        -s -name \\
        -fo ${sample_id}.original_cds.fa
    
    # Step 2: apply SNVs to generate mutated genome FASTA
    bcftools consensus \\
        -f ${reference_fasta} \\
        -o genome_mutated.fa \\
        ${filtered_vcf} > genome_mutated.fa

    # Step 3: extract mutated CDS from the mutated genome
    bedtools getfasta \\
        -fi genome_mutated.fa \\
        -bed \${frame_bed} \\
        -s -name \\
        -fo ${sample_id}.mutated_cds.fa

    # QC Check
    echo "Sample: ${sample_id}"
    echo "Original CDS count: \$(grep -c '^>' ${sample_id}.original_cds.fa)"
    echo "Mutated CDS count:  \$(grep -c '^>' ${sample_id}.mutated_cds.fa)"
    """
}

workflow {
    take:
        input_ch
    main:
        result = extract_sequences(input_ch)
    emit:
        original_cds_fa   = result.original_cds_fa
        genome_mutated_fa = result.genome_mutated_fa
        mutated_cds_fa    = result.mutated_cds_fa
}
