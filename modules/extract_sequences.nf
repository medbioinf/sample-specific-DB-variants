// Extracts coding sequences (CDS) from the reference genome and from a sample-specific genome containing protein-altering mutations.
// Requires 5 items:
//      reference genome: .fa + .fa.fai
//      BED file with CDS regions of affected transcripts (extract_cds_bed.nf output): .bed 
//      protein altering variants (protein_altering_bundle output): .vcf.gz + .vcf.gz.tbi

process extract_sequences {
    tag "$sample_id"
    publishDir "${params.outdir}/vcf/extract_sequences", mode: 'copy', overwrite: true, dynamic: true
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

    # Step 1: extract original CDS from reference genome
    bedtools getfasta \\
        -fi ${reference_fasta} \\
        -bed ${cds_affected_bed} \\
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
        -bed ${cds_affected_bed} \\
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