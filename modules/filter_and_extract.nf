// 1) Filter annotated VCF for protein-altering variants 
// 2) Extracts the transcripts IDs where the protein-altering variants occur
// Requires 3 items in the input tuple: 
//      sample_id
//      annotated variants .vcf.gz + .tbi 

process filter_and_extract {
    tag "$sample_id"
    publishDir "${params.outdir}/vcf/protein_altering", mode: 'copy', overwrite: true, dynamic: true
    cpus { params.threads ?: 8 }
    container (params.toolbox_image ?: 'biocontainers/bcftools:v1.20-1--deb_cv1')
    
    input:
    tuple val(sample_id), path(annotated_vcf), path(tbi_file)

    when:
    annotated_vcf.exists() && tbi_file.exists()

    output:
    tuple val(sample_id),
       path("${sample_id}.protein_altering.vcf.gz"),
       path("${sample_id}.protein_altering.vcf.gz.tbi"),
       path("${sample_id}.protein_altering.positions.tsv"),
       path("${sample_id}.protein_altering.bed"),
       path("${sample_id}.affected_transcripts.txt"),
       emit: protein_altering_bundle
    
    script:
    """
    echo "Filtering protein-altering variants for sample: ${sample_id}"

    # Step 1: Extract lines with protein-altering consequences from BCSQ
    # Keep: missense/frameshift/stop_gained/stop_lost/start_lost 
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/BCSQ\\n' ${annotated_vcf} \\
        | grep -E 'missense|frameshift|stop_gained|stop_lost|start_lost' \\
        | cut -f1,2 > ${sample_id}.protein_altering.positions.tsv

    # Step 2: Convert to BED format (0-based)
    awk '{print \$1"\\t"\$2-1"\\t"\$2}' ${sample_id}.protein_altering.positions.tsv > ${sample_id}.protein_altering.bed

    # Step 3: Filter annotated VCF based on the BED regions
    bcftools view -R ${sample_id}.protein_altering.bed ${annotated_vcf} -Oz -o ${sample_id}.protein_altering.vcf.gz

    # Step 4: Index filtered VCF
    bcftools index -f -t ${sample_id}.protein_altering.vcf.gz
    
    echo "Filtering completed for sample: ${sample_id}"

    # Step 5: Extract affected transcript IDs from BCSQ (unique ENST*)
    bcftools query -f '%INFO/BCSQ\\n' ${sample_id}.protein_altering.vcf.gz \\
        | tr ',' '\\n' \\
        | awk -F'|' '{
            tx=\"\"
            for (i=1; i<=NF; i++) {
                if (\$i ~ /^ENST/) { tx=\$i; break }
                if (\$i ~ /^transcript:/) {
                    split(\$i, parts, \":\")
                    if (parts[2] ~ /^ENST/) { tx=parts[2]; break }
                }
            }
            if (tx != \"\") print tx
        }' \\
        | sort -u > ${sample_id}.affected_transcripts.txt

    echo "Number of affected transcripts: \$(wc -l < ${sample_id}.affected_transcripts.txt)"
    """
}

// Workflow definition
workflow {
  take:
    input_ch  // tuple(val(sample_id), path(annotated.vcf.gz), path(annotated.vcf.gz.tbi))

  main:
    out_bundle = filter_and_extract(input_ch)

  emit:
    protein_altering_bundle = out_bundle.protein_altering_bundle
}
