process gatk_haplotypecaller {
  tag "$sample_id"
  publishDir "${params.outdir}/vcf/gatk", mode: 'copy', overwrite: true
  cpus { params.threads ?: 8 }
  memory '12 GB'
  container (params.toolbox_image ?: 'broadinstitute/gatk:4.5.0.0')
  // conda 'bioconda::gatk4=4.5.0.0'

  input:
    tuple val(sample_id), path(split_bam)      // from split_n_cigar (GATK path only)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)      // from ch_ref_for_gatk
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

  # GATK HaplotypeCaller -> RAW VCF (bgzipped) and index
  gatk HaplotypeCaller \
    -R ${ref_fasta} \
    -I ${split_bam} \
    -O ${sample_id}.raw.vcf.gz \
    --dont-use-soft-clipped-bases true \
    --standard-min-confidence-threshold-for-calling 20.0 \
    --create-output-variant-index true

  # Ensure tabix index exists (GATK should create it; make it explicit just in case)
  [ -f ${sample_id}.raw.vcf.gz.tbi ] || tabix -p vcf ${sample_id}.raw.vcf.gz
  
  # DP/ALT filter for GATK VCFs
  # - total depth: INFO/DP (site) or FORMAT/DP (sample)
  # - ALT reads:   FORMAT/AD[0:1]  (sample 0, first ALT)
  bcftools norm -m -any -Ou ${sample_id}.raw.vcf.gz \
  | bcftools view -i '((INFO/DP >= ${min_dp}) || (FORMAT/DP >= ${min_dp})) \
                    && (FORMAT/AD[0:1] >= ${min_alt})' \
    -Oz -o ${sample_id}.filtered.vcf.gz

  bcftools index -t ${sample_id}.filtered.vcf.gz
  """
}