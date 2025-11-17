
process sort_index_bam {
  tag "$sample_id"
  publishDir "${params.outdir}/04_bam", mode: 'copy', overwrite: true
  cpus (params.threads ?: 8)
  memory '8 GB'
  container (params.toolbox_image ?: 'biocontainers/samtools:v1.20-1-deb_cv1')
  // conda 'bioconda::samtools=1.20'

  input:
    tuple val(sample_id), path(aln_in)   // ${sample}.sam OR ${sample}.unsorted.bam

  output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam
    path "${sample_id}.sorted.bam.bai"

  script:
  """
  set -euo pipefail

  # samtools sort auto-detects SAM/BAM/CRAM input; default output is BAM
  samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${aln_in}

  # index
  samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
  """
}
