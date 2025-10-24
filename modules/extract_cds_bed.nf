// Isolates the CDS (coding regions) of only the transcripts affected by variants
// Requires 2 items as input:
//      Transcript annotation file gff3.gz + gff3.gz.tbi
//      Transcripts IDs where the protein-altering variants occur affected_transcripts.txt (extract_transcripts.nf output)

process extract_cds_bed {
    tag "$sample_id"
    publishDir "${params.outdir}/vcf/extract_cds_bed", mode: 'copy', overwrite: true, dynamic: true
    cpus { params.threads ?: 8 }
    container (params.toolbox_image ?: 'biocontainers/bcftools:v1.20-1--deb_cv1')

    input:
    tuple val(sample_id),
       path(positions_tsv),
       path(prot_alt_bed),
       path(filtered_vcf),
       path(filtered_vcf_tbi),
       path(affected_tx),
       path(gff3_file),
       path(gff3_tbi)

    when:
    positions_tsv.exists() && prot_alt_bed.exists() &&
    filtered_vcf.exists() && filtered_vcf_tbi.exists() &&
    affected_tx.exists() && gff3_file.exists() && gff3_tbi.exists()

    output:
    tuple val(sample_id),
       path("${sample_id}.protein_altering.positions.tsv"),
       path("${sample_id}.protein_altering.bed"),
       path("${sample_id}.protein_altering.vcf.gz"),
       path("${sample_id}.protein_altering.vcf.gz.tbi"),
       path("${sample_id}.affected_transcripts.txt"),
       path("${sample_id}.cds_affected.bed"),
       emit: cds_affected_bed

    script:
    """
    set -euo pipefail

    # Normalize input names so downstream stays consistent (avoid self-referential symlinks)
    for spec in \
      "${positions_tsv}:${sample_id}.protein_altering.positions.tsv" \
      "${prot_alt_bed}:${sample_id}.protein_altering.bed" \
      "${filtered_vcf}:${sample_id}.protein_altering.vcf.gz" \
      "${filtered_vcf_tbi}:${sample_id}.protein_altering.vcf.gz.tbi" \
      "${affected_tx}:${sample_id}.affected_transcripts.txt"
    do
      src="\${spec%%:*}"
      dest="\${spec##*:}"
      if [ "\${src}" != "\${dest}" ]; then
        ln -sf "\${src}" "\${dest}"
      fi
    done
    
    # Step 1: Extract CDS regions from GFF3 that match affected transcripts
    zcat -f "${gff3_file}" | grep -w "CDS" | \
    awk -v ids_file="${sample_id}.affected_transcripts.txt" '
      BEGIN{
        FS=OFS="\t"
        # load transcript IDs we want to keep
        while ((getline t < ids_file) > 0) keep[t]=1
      }
      {
        # extract transcript ID from Parent attribute (two common styles)
        id=""
        if (match(\$9, /Parent=transcript:([^;]+)/, m)) id=m[1]
        else if (match(\$9, /Parent=([^;]+)/, m))       id=m[1]

        # keep CDS only if its transcript is in the affected set
        if (id in keep) {
            s=\$4-1; if (s<0) s=0  # 0-based start
            print \$1, s, \$5, id, ".", \$7
        }
      }
    ' > "${sample_id}.cds_affected.bed"

    # Step 2: Print summary and preview for debugging
    echo "Total extracted CDS rows for ${sample_id}: \$(wc -l < ${sample_id}.cds_affected.bed)"
    echo "Preview of ${sample_id}.cds_affected.bed:"
    head -n 20 "${sample_id}.cds_affected.bed" || true
    """
}

// Workflow definition
workflow {
  take:
    input_ch  // tuple(sid, positions.tsv, bed, vcf.gz, tbi, affected_tx, gff3.gz, gff3.gz.tbi)

  main:
    res = extract_cds_bed(input_ch)

  emit:
    cds_affected_bed = res.cds_affected_bed
}
