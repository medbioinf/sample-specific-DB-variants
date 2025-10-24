// Annotate variants
// Requires a tuple that contains 7 items: 
//          sample_id
//          .vcf + vcf.tbi  
//          reference .fa +.fa.fai 
//          gff3 + gff3.tbi

process annotate_variants {
  tag "$sample_id"
  publishDir "${params.outdir}/vcf/annotated_variants", mode: 'copy'
  cpus { params.threads ?: 8 }
  container (params.toolbox_image ?: 'biocontainers/bcftools:v1.20-1--deb_cv1')

  input:
  tuple val(sample_id),
    path(snv_vcf),
    path(snv_tbi),
    path(reference_fasta),
    path(reference_fai),
    path(gff3_file),
    path(gff3_tbi)
  
  when:
  snv_vcf.exists() && snv_tbi.exists() &&
  reference_fasta.exists() && reference_fai.exists() &&
  gff3_file.exists() && gff3_tbi.exists()

  output:
  tuple val(sample_id),
    path("${sample_id}.annotated.vcf.gz"),
    path("${sample_id}.annotated.vcf.gz.tbi"),
    emit: annotated_vcf

  script:
  """
  set -euo pipefail
  bcftools csq -f ${reference_fasta} -g ${gff3_file} ${snv_vcf} -p a -Oz -o ${sample_id}.annotated.vcf.gz
  bcftools index -f -t ${sample_id}.annotated.vcf.gz
  """
}