// Identify A-to-I RNA editing events among protein-altering variants
// Requires:
//    sample_id
//    protein-altering VCF + index
//    BED/TSV of known editing sites (at least chrom + pos, optional end)

process match_rna_editing {
  tag "$sample_id"
  publishDir "${params.outdir}/06_vcf/03_rna_editing_matches", mode: 'copy', overwrite: true, dynamic: true
  cpus { params.threads ?: 4 }
  container (params.toolbox_image ?: 'ghcr.io/your-org/variant-caller-toolbox:latest')

  input:
  tuple val(sample_id),
        path(protein_vcf),
        path(protein_vcf_tbi),
        path(editing_catalog)

  when:
  protein_vcf.exists() && protein_vcf_tbi.exists() && editing_catalog.exists()

  output:
  tuple val(sample_id),
        path("${sample_id}.rna_editing_matches.txt"),
        emit: editing_matches

  script:
  """
  set -euo pipefail

  if [ ! -s "${editing_catalog}" ]; then
    echo "RNA editing catalog appears empty. Emitting placeholder matches for ${sample_id}."
    : > ${sample_id}.rna_editing_matches.txt
    exit 0
  fi

  python3 - <<'PY' "${protein_vcf}" "${editing_catalog}" "${sample_id}.rna_editing_matches.txt"
import gzip
import sys
from collections import defaultdict

vcf_path, catalog_path, out_path = sys.argv[1:]

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def load_catalog(path):
    sites = defaultdict(set)
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                if len(parts) >= 3:
                    start = int(parts[1])
                    pos = start + 1  # assume BED-style (0-based, half-open)
                else:
                    pos = int(parts[1])
            except ValueError:
                continue
            sites[chrom].add(pos)
    return sites

try:
    catalog_sites = load_catalog(catalog_path)
except Exception as exc:
    sys.stderr.write(f"[match_rna_editing] Failed to read catalog '{catalog_path}': {exc}\\n")
    sys.exit(1)

with open(out_path, "w") as out_handle:
    first = True
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as vcf_handle:
        for line in vcf_handle:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split("\\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            ref = parts[3].strip().upper()
            alts = [a.strip().upper() for a in parts[4].split(",") if a.strip()]
            chrom_sites = catalog_sites.get(chrom)
            if not chrom_sites or pos not in chrom_sites or len(ref) != 1:
                continue
            for alt in alts:
                if len(alt) != 1:
                    continue
                if (ref == "A" and alt == "G") or (ref == "T" and alt == "C"):
                    if first:
                        out_handle.write("CHROM\\tPOS\\tREF\\tALT\\n")
                        first = False
                    out_handle.write(f"{chrom}\\t{pos}\\t{ref}\\t{alt}\\n")

if first:
    with open(out_path, "w") as out_handle:
        out_handle.write("CHROM\\tPOS\\tREF\\tALT\\n")
PY

  """
}

workflow {
  take:
    input_ch   // tuple(val(sample_id), path(protein_vcf), path(protein_vcf_tbi), path(editing_catalog))

  main:
    matches = match_rna_editing(input_ch)

  emit:
    editing_matches = matches.editing_matches
}
