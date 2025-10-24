// Translates CDS into protein sequences 
// Writes a per-peptide annotation file
// Requires 2 items:
//      Sample ID: sample_id
//      FASTA with CDS entries of variants-affected transcripts: .mutated_cds.fa (extract_sequences.nf output)
//      Filtered variants that affect CDS: .filtered.vcf + .vcf.tbi (filter_variants.nf output)

process translate_proteins {
    tag "$sample_id"
    publishDir "${params.outdir}/vcf/translate_proteins", mode: 'copy'
    cpus { params.threads ?: 8 }
    container (params.toolbox_image ?: 'ghcr.io/your-org/variant-caller-toolbox:latest')

    input:
    tuple val(sample_id),
      path(mutated_cds_fa),
      path(filtered_vcf),
      path(filtered_vcf_tbi)

    when:
    mutated_cds_fa.exists() && filtered_vcf.exists() && filtered_vcf_tbi.exists()

    output:
    path("${sample_id}.merged_cds.fa"),         emit: merged_cds_fa
    path("${sample_id}.mutated_proteins.fa"),   emit: mutated_proteins_fa
    path("${sample_id}.mutated_proteins.annot.txt"), emit: mutated_proteins_annot

    script:
    """
    set -euo pipefail

    # Step 1: merge exon chunks per transcript (headers may be '>ID' or '>ID:...')
    awk '
      BEGIN { id=""; n=0 }
      /^>/ {
        hdr=substr(\$0,2)
        split(hdr,parts,":")
        t=parts[1]
        if(!(t in seen)){
          seen[t]=1
          order[++n]=t
        }
        id=t
        next
      }
      {
        seq[id]=seq[id]\$0
      }
      END {
        for(i=1;i<=n;i++){
          t=order[i]
          print ">" t
          print seq[t]
        }
      }
    ' ${mutated_cds_fa} > ${sample_id}.merged_cds.fa
    
    # check
    if [ ! -s ${sample_id}.merged_cds.fa ]; then
      echo "Error: ${sample_id}.merged_cds.fa is empty or missing" >&2
      exit 2
    fi

    # Step 2: translate merged CDS
    # genetic code (NCBI table) can be overridden via --genetic_code (default 1)
    GC=${params.genetic_code ?: 1}
    transeq -sequence ${sample_id}.merged_cds.fa -outseq ${sample_id}.mutated_proteins.fa -frame 1 -table \$GC

    echo "Proteins for ${sample_id}: \$(grep -c '^>' ${sample_id}.mutated_proteins.fa || true)"

    # Step 3: annotate peptides with variant info
    # collect transcript IDs present in the merged CDS
    grep '^>' ${sample_id}.merged_cds.fa | sed 's/^>//; s/:.*//' > tx.list
    export CALLER="${params.caller}"
    export SAMPLE_ID="${sample_id}"

    # extract BCSQ once (CHROM,POS,BCSQ); BCSQ entries are comma-separated; fields split by '|'
    bcftools query -f '%CHROM\\t%POS\\t%INFO/BCSQ\\n' ${filtered_vcf} > bcsq.tsv

    python3 - <<'PY'
import csv
import collections
import os

aa_map = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "U": "Sec", "O": "Pyl"
}

def expand(code: str) -> str:
    if len(code) == 1:
        return aa_map.get(code, code)
    return ''.join(aa_map.get(ch, ch) for ch in code)

def parse_change(field: str):
    if ">" not in field:
        return None
    left, right = field.split(">", 1)
    pos = ''.join(ch for ch in left if ch.isdigit())
    ref = ''.join(ch for ch in left if ch.isalpha() or ch == "*")
    alt = ''.join(ch for ch in right if ch.isalpha() or ch == "*")
    if not pos or not ref or not alt:
        return None
    return pos, ref, alt

ann = collections.OrderedDict()
caller = os.environ.get("CALLER", "NA")
sample = os.environ.get("SAMPLE_ID", "sample")

with open("bcsq.tsv") as handle:
    reader = csv.reader(handle, delimiter="\\t")
    for row in reader:
        if len(row) < 3 or not row[2]:
            continue
        for entry in row[2].split(","):
            fields = entry.split("|")
            tx_id = ""
            for field in fields:
                if field.startswith("ENST"):
                    tx_id = field
                    break
                if field.startswith("transcript:"):
                    _, _, maybe = field.partition(":")
                    if maybe.startswith("ENST"):
                        tx_id = maybe
                        break
            if not tx_id:
                continue
            aa_change = None
            for field in fields:
                parsed = parse_change(field)
                if not parsed:
                    continue
                pos, ref, alt = parsed
                aa_change = f"p.{expand(ref)}{pos}{expand(alt)}"
                break
            if not aa_change:
                continue
            bucket = ann.setdefault(tx_id, [])
            if aa_change not in bucket:
                bucket.append(aa_change)

with open("tx2aa.tsv", "w") as out_handle:
    for tx_id, changes in ann.items():
        out_handle.write(f"{tx_id}\\t{';'.join(changes)}\\n")

with open("tx.list") as tx_handle, open(f"{sample}.mutated_proteins.annot.txt", "w") as annot_handle:
    for raw in tx_handle:
        tx = raw.strip()
        if not tx:
            continue
        variants = ann.get(tx, [])
        aa_str = ';'.join(variants) if variants else "NA"
        annot_handle.write(f">{tx}\\n")
        annot_handle.write(f"caller: {caller}\\n")
        annot_handle.write(f"variants: {aa_str}\\n")
PY

    echo "Proteins: \$(grep -c '^>' ${sample_id}.mutated_proteins.fa || true)"
    echo "Annotations written to ${sample_id}.mutated_proteins.annot.txt"
    """
}

// Workflow definition
workflow {
  take:
    input_ch   // tuple(val(sample_id), path(mutated_cds_fa), path(filtered_vcf), path(filtered_vcf_tbi))

  main:
    result = translate_proteins(input_ch)

  emit:
    mutated_proteins_fa    = result.mutated_proteins_fa
    mutated_proteins_annot = result.mutated_proteins_annot
}
