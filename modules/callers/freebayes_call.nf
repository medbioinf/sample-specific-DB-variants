process freebayes_call {
  tag "$sample_id"
  publishDir "${params.outdir}/05_var_ident/freebayes", mode: 'copy', overwrite: true
  cpus { params.threads ?: 8 }
  container (params.toolbox_image ?: 'ghcr.io/your-org/variant-caller-toolbox:latest')

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

  # freebayes -> RAW VCF (bgzipped) and index
  freebayes \
    -f ${ref_fasta} \
    --min-base-quality 20 \
    --min-mapping-quality 20 \
    --ploidy 2 \
    ${bam_for_call} \
  | bgzip -c -@ ${task.cpus} > ${sample_id}.raw.vcf.gz

  tabix -p vcf ${sample_id}.raw.vcf.gz

  # DP/ALT filter -> filtered VCF and index
  # INFO/DP - total depth at the site, both ref and alt alleles
  # FORMAT/AO or FORMAT/SAF + FORMATE/SAR - depth of alt alleles on forward/reverse strand
  bcftools view -i 'INFO/DP >= ${min_dp} && ( \
                    FORMAT/AO[0:0] >= ${min_alt} || \
                    (INFO/SAF + INFO/SAR) >= ${min_alt} )' \
  ${sample_id}.raw.vcf.gz -Oz -o ${sample_id}.filtered.vcf.gz

  bcftools index -t ${sample_id}.filtered.vcf.gz

  echo "Variants calling step is completed"
  """
}
