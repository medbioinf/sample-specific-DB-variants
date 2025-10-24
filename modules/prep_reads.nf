process prep_reads {
  tag "$sample_id"
  publishDir "${params.outdir}/qc/fastp", mode: 'copy', overwrite: true
  cpus (params.threads ?: 8)
  memory '4 GB'
  container (params.toolbox_image ?: 'quay.io/biocontainers/fastp:0.23.4--h5f7400f_0')
  // conda 'bioconda::fastp=0.23.4'

  input:
    // Works with either a list [R1,R2] or a single file
    tuple val(sample_id), path(reads)

  output:
    // Emit one or two trimmed FASTQs; downstream modules just take the tuple
    tuple val(sample_id), path("${sample_id}_R*.trim.fastq.gz"), emit: trimmed
    path "${sample_id}.fastp.html"
    path "${sample_id}.fastp.json"

  script:
  // Normalize to array
  def rs = (reads instanceof List) ? reads : [ reads ]
  def r1 = rs.find { it.name =~ /R1|_1\./ } ?: rs[0]
  def r2 = rs.find { it.name =~ /R2|_2\./ }

  """
  set -euo pipefail

  # --- Hardcoded-but-sane defaults for RNA-seq ---
  CUT_QUAL=20             # mean quality for sliding trimming
  MIN_LEN=35              # drop very short reads after trim
  THREADS=${task.cpus}

  if [ -n "${r2 ?: ''}" ]; then
    # Paired-end
    fastp \
      -i ${r1} -I ${r2} \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --cut_front --cut_tail --cut_mean_quality \$CUT_QUAL \
      --length_required \$MIN_LEN \
      --thread \$THREADS \
      --html ${sample_id}.fastp.html \
      --json ${sample_id}.fastp.json \
      --report_title "fastp: ${sample_id}"
  else
    # Single-end
    fastp \
      -i ${r1} \
      -o ${sample_id}_SE.trim.fastq.gz \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --cut_front --cut_tail --cut_mean_quality \$CUT_QUAL \
      --length_required \$MIN_LEN \
      --thread \$THREADS \
      --html ${sample_id}.fastp.html \
      --json ${sample_id}.fastp.json \
      --report_title "fastp: ${sample_id}"
  fi
  """
}
