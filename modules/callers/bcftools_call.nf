process bcftools_call {
  tag "$sample_id"
  publishDir "${params.outdir}/vcf/bcftools", mode: 'copy', overwrite: true
  cpus { params.threads ?: 8 }
  container (params.toolbox_image ?: 'biocontainers/bcftools:v1.20-1--deb_cv1')
  // conda 'bioconda::bcftools=1.20'

  input:
    tuple val(sample_id), path(bam_for_call)   // markdup BAM (not split)
    path ref_fasta
    val  min_dp
    val  min_alt

  output:
    tuple val(sample_id), path("${sample_id}.raw.vcf.gz"),       emit: raw_vcf_gz
    path "${sample_id}.raw.vcf.gz.tbi",                          emit: raw_tbi
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"),  emit: vcf_gz
    path "${sample_id}.filtered.vcf.gz.tbi",                     emit: tbi

  script:
  """
  set -euo pipefail

  # mpileup -> call -> RAW VCF and index
  bcftools mpileup \
    -f ${ref_fasta} \
    -q 20 -Q 20 \
    -Ou \
    -a AD,DP,SP,ADF,ADR,DP4 \
    ${bam_for_call} \
  | bcftools call \
    -m -v \
    --ploidy 2 \
    -Oz -o ${sample_id}.raw.vcf.gz
  
  bcftools index -t ${sample_id}.raw.vcf.gz

  # DP/ALT filter (bcftools annotations) -> FILTERed VCF and index
  # INFO/DP - total depth at the site, both ref and alt alleles
  # FORMAT/AD or FORMAT/ADF + FORMATE/ADR - depth of alt alleles on forward/reverse strand
  bcftools view -i 'INFO/DP >= ${min_dp} && ( \
                    FORMAT/AD[0:1] >= ${min_alt} || \
                    (FORMAT/ADF[0:1] + FORMAT/ADR[0:1]) >= ${min_alt} || \
                    (INFO/DP4[2] + INFO/DP4[3]) >= ${min_alt} )' \
  ${sample_id}.raw.vcf.gz -Oz -o ${sample_id}.filtered.vcf.gz

  bcftools index -t ${sample_id}.filtered.vcf.gz
  """
}