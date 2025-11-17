// Translates CDS into protein sequences 
// Writes a per-peptide annotation file
// Requires 2 items:
//      Sample ID: sample_id
//      FASTA with CDS entries of variants-affected transcripts: .mutated_cds.fa (extract_sequences.nf output)
//      Filtered variants that affect CDS: .filtered.vcf + .vcf.tbi (filter_variants.nf output)

process translate_proteins {
    tag "$sample_id"
    publishDir "${params.outdir}/07_translate_proteins", mode: 'copy'
    cpus { params.threads ?: 8 }
    container (params.toolbox_image ?: 'ghcr.io/your-org/variant-caller-toolbox:latest')

    input:
    tuple val(sample_id),
      path(mutated_cds_fa),
      path(filtered_vcf),
      path(filtered_vcf_tbi),
      path(rna_editing_matches)

    when:
    mutated_cds_fa.exists() && filtered_vcf.exists() && filtered_vcf_tbi.exists()

    output:
    path("${sample_id}.merged_cds.fa"),         emit: merged_cds_fa
    path("${sample_id}.mutated_proteins.fa"),   emit: mutated_proteins_fa
    path("${sample_id}.mutated_proteins.annot.txt"), emit: mutated_proteins_annot

    script:
    """
    set -euo pipefail

    if [ ! -s ${mutated_cds_fa} ]; then
      echo "No CDS sequences to translate for sample: ${sample_id}"
      : > ${sample_id}.merged_cds.fa
      : > ${sample_id}.mutated_proteins.fa
      : > ${sample_id}.mutated_proteins.annot.txt
      exit 0
    fi

    # Step 1: merge exon chunks per transcript (headers may be '>ID' or '>ID:...')
    awk '
      BEGIN { id=""; n=0 }
      /^>/ {
        hdr=substr(\$0,2)
        split(hdr,parts,":")
        t=parts[1]
        if(!(t in seen)){
          seen[t]=1
          order[++n]=t
        }
        id=t
        next
      }
      {
        seq[id]=seq[id]\$0
      }
      END {
        for(i=1;i<=n;i++){
          t=order[i]
          print ">" t
          print seq[t]
        }
      }
    ' ${mutated_cds_fa} > ${sample_id}.merged_cds.fa
    
    # check
    if [ ! -s ${sample_id}.merged_cds.fa ]; then
      echo "Error: ${sample_id}.merged_cds.fa is empty or missing" >&2
      exit 2
    fi

    # Step 2: translate merged CDS
    # genetic code (NCBI table) can be overridden via --genetic_code (default 1)
    GC=${params.genetic_code ?: 1}
    transeq -sequence ${sample_id}.merged_cds.fa -outseq ${sample_id}.mutated_proteins.fa -frame 1 -table \$GC

    echo "Proteins for ${sample_id}: \$(grep -c '^>' ${sample_id}.mutated_proteins.fa || true)"

    # Step 3: annotate peptides with variant info
    # collect transcript IDs present in the merged CDS
    grep '^>' ${sample_id}.merged_cds.fa | sed 's/^>//; s/:.*//' > tx.list
    export CALLER="${params.caller}"
    export SAMPLE_ID="${sample_id}"

    # extract BCSQ once (CHROM,POS,BCSQ); BCSQ entries are comma-separated; fields split by '|'
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/BCSQ\\n' ${filtered_vcf} > bcsq.tsv

    if [ -n "\${rna_editing_matches:-}" ] && [ -s "\${rna_editing_matches:-}" ]; then
      export RNA_EDITING_MATCHES="${rna_editing_matches}"
    else
      unset RNA_EDITING_MATCHES || true
    fi

    python3 - <<'PY'
import collections
import csv
import os
import re

caller = os.environ.get("CALLER", "NA")
sample = os.environ.get("SAMPLE_ID", "sample")
editing_file = os.environ.get("RNA_EDITING_MATCHES")

editing_hits = set()
if editing_file and os.path.exists(editing_file):
    with open(editing_file) as matches:
        for line in matches:
            line = line.strip()
            if not line or line.startswith("CHROM"):
                continue
            parts = line.split("\\t")
            if len(parts) < 4:
                continue
            chrom, pos, _ref, alt = parts[:4]
            try:
                editing_hits.add((chrom, int(pos), alt.upper()))
            except ValueError:
                continue

ann = collections.OrderedDict()

with open("bcsq.tsv") as handle:
    reader = csv.reader(handle, delimiter="\\t")
    for row in reader:
        if len(row) < 5:
            continue
        chrom, pos_str, ref, alt_field, info = row
        if not info:
            continue
        try:
            pos = int(pos_str)
        except ValueError:
            continue
        entries = [entry for entry in info.split(",") if entry]
        for entry in entries:
            fields = entry.split("|")
            if len(fields) < 4:
                continue
            tx_id = fields[2] if len(fields) > 2 else ""
            if not tx_id.startswith("ENST"):
                continue
            consequence = fields[0]
            aa_change = fields[5] if len(fields) > 5 else ""
            dna_change = fields[6] if len(fields) > 6 else ""
            is_utr = "utr" in consequence.lower()
            if not aa_change or aa_change in {".", "-"}:
                if is_utr:
                    aa_change = "variant hits the UTR - no amino acid change"
                else:
                    continue
            alt_nt = None
            if dna_change and ">" in dna_change:
                match = re.search(r'>\\s*([ACGT]+)', dna_change)
                if match:
                    alt_nt = match.group(1).upper()
            is_editing = False
            if alt_nt:
                if (chrom, pos, alt_nt) in editing_hits:
                    is_editing = True
            bucket = ann.setdefault(tx_id, [])
            record = (aa_change, is_editing)
            if record not in bucket:
                bucket.append(record)

with open("tx2aa.tsv", "w") as out_handle:
    for tx_id, changes in ann.items():
        labels = []
        for aa_change, is_editing in changes:
            label = aa_change
            if is_editing:
                label = f"{label}(A->I editing)"
            labels.append(label)
        out_handle.write(f"{tx_id}\\t{';'.join(labels)}\\n")

with open("tx.list") as tx_handle, open(f"{sample}.mutated_proteins.annot.txt", "w") as annot_handle:
    for raw in tx_handle:
        tx = raw.strip()
        if not tx:
            continue
        variants = ann.get(tx, [])
        aa_str = ';'.join(
            f"{aa_change}{'(A->I editing)' if is_editing else ''}"
            for aa_change, is_editing in variants
        ) if variants else "NA"
        annot_handle.write(f">{tx}\\n")
        annot_handle.write(f"caller: {caller}\\n")
        annot_handle.write(f"variants: {aa_str}\\n")
PY

    echo "Proteins: \$(grep -c '^>' ${sample_id}.mutated_proteins.fa || true)"
    echo "Annotations written to ${sample_id}.mutated_proteins.annot.txt"
    echo "CDSs are translated to proteins"
    """
}

// Workflow definition
workflow {
  take:
    input_ch   // tuple(val(sample_id), path(mutated_cds_fa), path(filtered_vcf), path(filtered_vcf_tbi), path(rna_editing_matches?))

  main:
    result = translate_proteins(input_ch)

  emit:
    mutated_proteins_fa    = result.mutated_proteins_fa
    mutated_proteins_annot = result.mutated_proteins_annot
}
