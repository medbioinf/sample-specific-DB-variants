process fixmate {
  tag "$sample_id"
  publishDir "${params.outdir}/04_bam", mode: 'copy', overwrite: true
  cpus (params.fixmate_cpus ?: params.threads ?: 4)
  memory { params.fixmate_memory ?: '24 GB' }
  container (params.toolbox_image ?: 'biocontainers/samtools:v1.20-1-deb_cv1')

  input:
    tuple val(sample_id), path(sorted_bam)  // coordinate-sorted BAM from sort_index_bam

  output:
    tuple val(sample_id), path("${sample_id}.fixmate.cs.bam"), emit: bam  // coordinate-sorted with fixmate tags

  def sortMem = params.fixmate_sort_mem ?: '1G'

  script:
  """
  set -euo pipefail

  # name-sort, required by fixmate (limit per-thread memory to avoid tmp overflow)
  samtools sort -n -@ ${task.cpus} -m ${sortMem} \
    -T ${sample_id}.name_sorted.tmp \
    -o ${sample_id}.name_sorted.bam ${sorted_bam}

  # add mate tags
  samtools fixmate -m -@ ${task.cpus} \
    ${sample_id}.name_sorted.bam ${sample_id}.fixmate.ns.bam

  # re-sort to coordinate for markdup
  samtools sort -@ ${task.cpus} -m ${sortMem} \
    -T ${sample_id}.fixmate.cs.tmp \
    -o ${sample_id}.fixmate.cs.bam ${sample_id}.fixmate.ns.bam
  """
}
