process extend_protein_db_proc {
  tag "$sample_id"
  publishDir "${params.outdir}/08_final_db", mode: 'copy', overwrite: true, dynamic: true
  cpus { params.threads ?: 2 }
  container (params.toolbox_image ?: 'ghcr.io/your-org/variant-caller-toolbox:latest')

  input:
  tuple val(sample_id),
        path(mutated_proteins_fa),
        path(annotation_file),
        path(canonical_fasta),
        path(decoy_fasta)

  output:
  tuple val(sample_id),
        path("${sample_id}.extended_proteins.fa"),
        path("${sample_id}.mutated_proteins.annot.txt"),
        emit: extended_protein_db

  when:
    mutated_proteins_fa.exists() &&
    annotation_file.exists() &&
    canonical_fasta.exists() &&
    decoy_fasta.exists()

  script:
  """
  set -euo pipefail

  clean_mut=${sample_id}.mutated_proteins.cleaned.fa
  if [ -s "${mutated_proteins_fa}" ]; then
    awk 'NR==1 && \$0 !~ /^>/ { print ">"\$0; next } { print }' ${mutated_proteins_fa} > \${clean_mut}
  else
    : > \${clean_mut}
  fi

  python3 - <<'PY' "${sample_id}.extended_proteins.fa" "\${clean_mut}" "${canonical_fasta}" "${decoy_fasta}"
import sys
from pathlib import Path

out_path = Path(sys.argv[1])
fasta_files = [Path(p) for p in sys.argv[2:]]

with out_path.open("w") as out_handle:
    first_entry = True
    for fasta_path in fasta_files:
        if not fasta_path.exists() or fasta_path.stat().st_size == 0:
            continue
        with fasta_path.open() as fasta_handle:
            for line in fasta_handle:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith(">"):
                    if not first_entry:
                        out_handle.write("\\n")
                    out_handle.write(line + "\\n")
                    first_entry = False
                else:
                    seq = line.replace("\\r", "").replace(" ", "").replace("\\t", "").replace("*", "")
                    out_handle.write(seq + "\\n")
PY

  if [ "${annotation_file}" != "${sample_id}.mutated_proteins.annot.txt" ]; then
    cp ${annotation_file} ${sample_id}.mutated_proteins.annot.txt
  else
    # annotation already staged with desired name; ensure timestamp updated
    touch ${sample_id}.mutated_proteins.annot.txt
  fi

  echo "Extended protein DB for ${sample_id} written to ${sample_id}.extended_proteins.fa"
  echo "Canincal protein database is extended with mutated proteins"
  """
}

workflow final_protein_db {
  take:
    input_ch

  main:
    extend_protein_db_proc(input_ch)

  emit:
    extended_protein_db = extend_protein_db_proc.out
}
