process split_n_cigar {
  tag "$sample_id"
  publishDir "${params.outdir}/bam", mode: 'copy'
  cpus (params.threads ?: 8)
  container (params.toolbox_image ?: 'variant-caller-toolbox:latest')
   // conda 'bioconda::gatk4=4.5.0.0 bioconda::samtools=1.20'

  input:
  tuple val(sample_id), path(markdup_bam)                // from markdup
  tuple path(ref_fasta), path(ref_fai), path(ref_dict)   // from ch_ref_for_gatk

  output:
  tuple val(sample_id), path("${sample_id}.split.bam"), emit: split_bam

  script:
  """
  set -euo pipefail

  # ensure input BAM is indexed (no-op if already indexed)
  [ -f ${markdup_bam}.bai ] || samtools index -@ ${task.cpus} ${markdup_bam}

  gatk SplitNCigarReads \
    -R ${ref_fasta} \
    -I ${markdup_bam} \
    -O ${sample_id}.split.bam

  samtools index -@ ${task.cpus} ${sample_id}.split.bam
  """
}