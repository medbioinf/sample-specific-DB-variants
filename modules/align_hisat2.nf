process align_hisat2 {
  tag "$sample_id"
  publishDir "${params.outdir}/alignment", mode: 'copy'
  cpus { params.threads ?: 8 }
  memory '12 GB'
  container "${params.toolbox_image}"

  input:
  tuple val(sample_id), path(trimmed_fastq)     // 1 or 2 FASTQs (SE or PE)
  path index_files

  output:
  tuple val(sample_id), path("${sample_id}.sam")

  script:
  """
  set -euo pipefail

  # Stage index (same as you already do)
  if [ -d "${index_files}" ]; then
    cp ${index_files}/*.ht2 .
  elif [[ "${index_files}" == *.tar.gz ]]; then
    tar -xzf ${index_files}
  else
    echo "ERROR: Unsupported index input: ${index_files}" >&2
    exit 3
  fi

  shopt -s nullglob
  arr=( *.1.ht2 )
  if [ \${#arr[@]} -eq 0 ]; then
    echo "ERROR: No *.1.ht2 found in work dir" >&2
    exit 1
  fi
  first="\${arr[0]}"
  idx_prefix="\${first%.1.ht2}"

  # Read groups
  RGID='${sample_id}'
  RGSM='${sample_id}'
  RGLB='lib1'
  RGPL='ILLUMINA'
  RGPU='${sample_id}.L001'
  RG_OPTS="--rg-id \$RGID --rg SM:\$RGSM --rg LB:\$RGLB --rg PL:\$RGPL --rg PU:\$RGPU"

  # Determine SE vs PE
  # Nextflow passes 'path(trimmed_fastq)' as an array when the glob matches 2 files.
  mapfile -t fq < <(printf "%s\n" ${trimmed_fastq})

  if [ \${#fq[@]} -eq 2 ]; then
    # PAIRED-END
    hisat2 -p ${task.cpus} -x "\$idx_prefix" \
      -1 "\${fq[0]}" -2 "\${fq[1]}" \
      --dta --no-unal \$RG_OPTS \
      -S ${sample_id}.sam
  else
    # SINGLE-END
    hisat2 -p ${task.cpus} -x "\$idx_prefix" \
      -U "\${fq[0]}" \
      --dta --no-unal \$RG_OPTS \
      -S ${sample_id}.sam
  fi
  """
}
