process markdup {
  tag "$sample_id"
  publishDir "${params.outdir}/04_bam", mode: 'copy', overwrite: true
  cpus (params.markdup_threads ?: params.threads ?: 4)
  memory '16 GB'
  container (params.toolbox_image ?: 'biocontainers/samtools:v1.20-1-deb_cv1')
  // conda 'bioconda::samtools=1.20'

  input:
    // coordinate-sorted BAM from sort_index_bam
    tuple val(sample_id), path(sorted_bam)

  output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), emit: bam
    path "${sample_id}.markdup.bam.bai"
    path "${sample_id}.markdup.metrics.txt"

  script:
  // mark only (do not remove); flip to '-r' if you need removal
  def rm_flag = (params.remove_dups as boolean) ? '-r' : ''

  """
  set -euo pipefail

  # samtools markdup expects coordinate-sorted input
  samtools markdup -@ ${task.cpus} ${rm_flag} \
                   -l 150 \
                   -f ${sample_id}.markdup.metrics.txt \
                   ${sorted_bam} ${sample_id}.markdup.bam

  samtools index -@ ${task.cpus} ${sample_id}.markdup.bam
  """
}
