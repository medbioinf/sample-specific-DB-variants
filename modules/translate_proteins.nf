// Translates CDS into protein sequences 
// Writes a per-peptide annotation file
// Requires 2 items:
//      Sample ID: sample_id
//      FASTA with CDS entries of variants-affected transcripts: .mutated_cds.fa (extract_sequences.nf output)
//      Filtered variants that affect CDS: .filtered.vcf + .vcf.tbi (filter_variants.nf output)

process translate_proteins {
    tag "$sample_id"
    publishDir "${params.outdir}/07_translate_proteins", mode: 'copy'
    cpus { params.threads ?: 8 }
    container (params.toolbox_image ?: 'ghcr.io/your-org/variant-caller-toolbox:latest')

    input:
    tuple val(sample_id),
      path(mutated_cds_fa),
      path(original_cds_fa),
      path(filtered_vcf),
      path(filtered_vcf_tbi),
      path(rna_editing_matches)

    when:
    mutated_cds_fa.exists() && original_cds_fa.exists() && filtered_vcf.exists() && filtered_vcf_tbi.exists()

    output:
    path("${sample_id}.merged_cds.fa"),         emit: merged_cds_fa
    path("${sample_id}.mutated_proteins.fa"),   emit: mutated_proteins_fa
    path("${sample_id}.mutated_proteins.annot.txt"), emit: mutated_proteins_annot

    script:
    """
    set -euo pipefail

    if [ ! -s ${mutated_cds_fa} ]; then
      echo "No CDS sequences to translate for sample: ${sample_id}"
      : > ${sample_id}.merged_cds.fa
      : > ${sample_id}.mutated_proteins.fa
      : > ${sample_id}.mutated_proteins.annot.txt
      exit 0
    fi

    merge_cds() {
      local infile=\$1
      local outfile=\$2
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
      ' "\${infile}" > "\${outfile}"
    }

    merge_cds ${mutated_cds_fa} ${sample_id}.merged_cds.fa
    merge_cds ${original_cds_fa} ${sample_id}.merged_original_cds.fa
    
    # check
    if [ ! -s ${sample_id}.merged_cds.fa ]; then
      echo "Error: ${sample_id}.merged_cds.fa is empty or missing" >&2
      exit 2
    fi

    # Step 2: translate merged CDS (mutated + original) to keep frame alignment
    GC=${params.genetic_code ?: 1}
    transeq -sequence ${sample_id}.merged_cds.fa -outseq ${sample_id}.mutated_proteins.fa -frame 1 -table \$GC
    transeq -sequence ${sample_id}.merged_original_cds.fa -outseq ${sample_id}.original_proteins.fa -frame 1 -table \$GC

    echo "Proteins for ${sample_id}: \$(grep -c '^>' ${sample_id}.mutated_proteins.fa || true)"

    # Step 3: annotate peptides with variant info
    export CALLER="${params.caller}"
    export SAMPLE_ID="${sample_id}"

    # extract BCSQ once (CHROM,POS,BCSQ); BCSQ entries are comma-separated; fields split by '|'
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/BCSQ\\n' ${filtered_vcf} > bcsq.tsv

    if [ -n "\${rna_editing_matches:-}" ] && [ -s "\${rna_editing_matches:-}" ]; then
      export RNA_EDITING_MATCHES="${rna_editing_matches}"
    else
      unset RNA_EDITING_MATCHES || true
    fi

    python3 - <<'PY'
import csv
import os
from collections import OrderedDict, defaultdict

caller = os.environ.get("CALLER", "NA")
sample = os.environ.get("SAMPLE_ID", "sample")
editing_file = os.environ.get("RNA_EDITING_MATCHES")


def normalize_id(raw):
    # Strip any trailing _<number> suffixes (one or more) so ENST00000372805_1_1 maps back to ENST00000372805.
    base = raw.split(":")[0]
    while True:
        parts = base.rsplit("_", 1)
        if len(parts) == 2 and parts[1].isdigit():
            base = parts[0]
            continue
        break
    return base

def read_fasta(path):
    order = []
    seqs = {}
    hdr_to_base = {}
    cur = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                raw_hdr = line[1:]
                cur = raw_hdr
                base = normalize_id(raw_hdr)
                hdr_to_base[cur] = base
                if cur not in order:
                    order.append(cur)
                if cur not in seqs:
                    seqs[cur] = []
            else:
                if cur is not None:
                    seqs.setdefault(cur, []).append(line)

    # Collapse duplicates by choosing the longest non-empty sequence for each id
    collapsed = {}
    for k, v in seqs.items():
        seq = "".join(v)
        if k not in collapsed or len(seq) > len(collapsed[k]):
            collapsed[k] = seq

    return order, collapsed, hdr_to_base

def write_clean_fasta(path, order, seqs, header_map=None):
    with open(path, "w") as handle:
        first = True
        for tx in order:
            seq = seqs.get(tx, "")
            if not seq:
                continue
            clean = seq.replace("*", "")
            if not first:
                handle.write("\\n")
            name = header_map.get(tx, tx) if header_map else tx
            handle.write(f">{name}\\n{clean}\\n")
            first = False

def compute_diffs(ref_seq, alt_seq):
    # Derive amino acid differences between reference and mutated proteins
    # Compare up to the shorter length; if lengths differ, append a single len change.
    diffs = []
    max_len = min(len(ref_seq), len(alt_seq))
    for i in range(max_len):
        a = ref_seq[i]
        b = alt_seq[i]
        if a != b:
            diffs.append(f"{i+1}{a}>{b}")
    if len(ref_seq) != len(alt_seq):
        diffs.append(f"len{len(ref_seq)}->{len(alt_seq)}")
    return diffs


editing_hits = set()
if editing_file and os.path.exists(editing_file):
    with open(editing_file) as matches:
        for line in matches:
            line = line.strip()
            if not line or line.startswith("CHROM"):
                continue
            parts = line.split("\\t")
            if len(parts) < 4:
                continue
            chrom, pos, _ref, alt = parts[:4]
            try:
                editing_hits.add((chrom, int(pos), alt.upper()))
            except ValueError:
                continue

# collect per-transcript variant details
variant_map = OrderedDict()
with open("bcsq.tsv") as handle:
    reader = csv.reader(handle, delimiter="\\t")
    for row in reader:
        if len(row) < 5:
            continue
        chrom, pos_str, ref, alt_field, info = row
        if not info:
            continue
        try:
            pos = int(pos_str)
        except ValueError:
            continue

        alts = [a for a in alt_field.split(",") if a]
        entries = [entry for entry in info.split(",") if entry]
        seen_keys = set()
        for idx, entry in enumerate(entries):
            fields = entry.split("|")
            if len(fields) < 3:
                continue
            tx_id = fields[2] if len(fields) > 2 else ""
            if not tx_id.startswith("ENST"):
                # try alternative transcript encoding
                for field in fields:
                    if field.startswith("transcript:") and field.split(":", 1)[1].startswith("ENST"):
                        tx_id = field.split(":", 1)[1]
                        break
            if not tx_id:
                continue
            gene = fields[1] if len(fields) > 1 else ""
            aa_change = fields[5] if len(fields) > 5 else ""
            dna_change = fields[6] if len(fields) > 6 else ""
            consequence = fields[0] if len(fields) > 0 else ""
            alt = alts[idx] if idx < len(alts) else (alts[-1] if alts else "")
            key = (chrom, pos, alt, tx_id, aa_change, gene)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            recs = variant_map.setdefault(tx_id, [])
            recs.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene,
                "aa": aa_change,
                "dna": dna_change,
                "cons": consequence,
                "is_editing": (chrom, pos, alt.upper()) in editing_hits if alt else False,
            })

# assign variant indices per transcript
for tx_id, records in variant_map.items():
    unique = []
    seen = set()
    for rec in records:
        key = (rec["chrom"], rec["pos"], rec["alt"], rec["aa"], rec["gene"])
        if key in seen:
            continue
        seen.add(key)
        unique.append(rec)
    for i, rec in enumerate(unique, 1):
        rec["idx"] = i
    variant_map[tx_id] = unique

# read translated proteins (mutated and reference); keep raw (with stops) for position adjustment
tx_order, prot_alt_raw, hdr_base_alt = read_fasta(f"{sample}.mutated_proteins.fa")
extra_order, prot_ref_raw, hdr_base_ref = read_fasta(f"{sample}.original_proteins.fa")
# ensure all tx from ref are included
for tx in extra_order:
    if tx not in tx_order:
        tx_order.append(tx)

# write cleaned FASTAs (without stop codons) for downstream consistency
write_clean_fasta(f"{sample}.mutated_proteins.fa", tx_order, prot_alt_raw)
write_clean_fasta(f"{sample}.original_proteins.fa", extra_order or tx_order, prot_ref_raw)

# Build cleaned dictionaries for diffing
clean_alt = {k: v.replace("*", "") for k, v in prot_alt_raw.items()}
clean_ref = {k: v.replace("*", "") for k, v in prot_ref_raw.items()}

with open(f"{sample}.mutated_proteins.annot.txt", "w") as annot_handle:
    # enforce 1:1 mapping: only headers present in FASTA, one block per header
    for tx in tx_order:
        base_tx = hdr_base_alt.get(tx) or hdr_base_ref.get(tx) or tx.split("_")[0]
        variants = variant_map.get(base_tx, [])
        header = tx
        annot_handle.write(f">{header}\\n")

        ref_clean_seq = clean_ref.get(tx, "")
        alt_clean_seq = clean_alt.get(tx, "")
        ref_raw_seq = prot_ref_raw.get(tx, "")
        alt_raw_seq = prot_alt_raw.get(tx, "")
        diffs = compute_diffs(ref_clean_seq, alt_clean_seq) if ref_clean_seq or alt_clean_seq else []
        aa_change_str = "; ".join(diffs) if diffs else "NA"
        stop_removed = (
            aa_change_str == "NA"
            and (ref_raw_seq or alt_raw_seq)
            and (ref_raw_seq.replace("*", "") == ref_clean_seq)
            and (alt_raw_seq.replace("*", "") == alt_clean_seq)
            and (ref_raw_seq != alt_raw_seq)
            and ("*" in ref_raw_seq or "*" in alt_raw_seq)
        )

        # flag A->I editing if any variant marked as editing
        has_editing = any(v.get("is_editing") for v in variants)
        if has_editing and aa_change_str != "NA":
            aa_change_str = f"{aa_change_str} (A->I editing)"
        elif has_editing and aa_change_str == "NA":
            aa_change_str = "A->I editing"

        if variants:
            # fallback: if no diff detected, use BCSQ AA or DNA change so we never drop info
            if aa_change_str == "NA":
                if stop_removed:
                    aa_change_str = "no amino acid change after stop removal (stop codon dropped)"
                else:
                    aa_fallback = [v.get("aa") for v in variants if v.get("aa")]
                    if aa_fallback:
                        aa_change_str = "; ".join(aa_fallback)
                    else:
                        dna_fallback = [v.get("dna") for v in variants if v.get("dna")]
                        if dna_fallback:
                            aa_change_str = "; ".join(dna_fallback)
                        else:
                            # include consequence to explain why no AA change (e.g., UTR/synonymous)
                            cons = list({v.get("cons") or "no_protein_change" for v in variants})
                            aa_change_str = "; ".join(f"no amino acid change ({c})" for c in cons)
            idx_line = "; ".join(str(v["idx"]) for v in variants)
            locus_parts = []
            for v in variants:
                alt = v.get("alt") or "NA"
                ref = v.get("ref") or "NA"
                locus_parts.append(f"{v['chrom']}:{v['pos']} {ref}>{alt}")
            locus_line = "; ".join(locus_parts)
            genes = []
            for v in variants:
                genes.append(v.get("gene") or "NA")
            annot_handle.write(f"variant_index: {idx_line}\\n")
            annot_handle.write(f"variant_locus: {locus_line}\\n")
            annot_handle.write(f"amino_acid_change: {aa_change_str}\\n")
            annot_handle.write(f"gene: {'; '.join(genes)}\\n")
        else:
            # No transcript hit in BCSQ: keep sequences aligned and report why
            if aa_change_str == "NA":
                if stop_removed:
                    aa_change_str = "no amino acid change after stop removal (stop codon dropped)"
                elif ref_clean_seq or alt_clean_seq:
                    aa_change_str = "no amino acid change (no mapped variant)"
                else:
                    aa_change_str = "no amino acid change (no sequence available)"
            annot_handle.write("variant_index: NA (no transcript match in BCSQ)\\n")
            annot_handle.write("variant_locus: NA (no transcript match in BCSQ)\\n")
            annot_handle.write(f"amino_acid_change: {aa_change_str}\\n")
            annot_handle.write("gene: NA (no transcript match in BCSQ)\\n")
        annot_handle.write(f"caller: {caller}\\n")
PY

    echo "Proteins: \$(grep -c '^>' ${sample_id}.mutated_proteins.fa || true)"
    echo "Annotations written to ${sample_id}.mutated_proteins.annot.txt"
    echo "CDSs are translated to proteins"
    """
}

// Workflow definition
workflow {
  take:
    input_ch   // tuple(val(sample_id), path(mutated_cds_fa), path(filtered_vcf), path(filtered_vcf_tbi), path(rna_editing_matches?))

  main:
    result = translate_proteins(input_ch)

  emit:
    mutated_proteins_fa    = result.mutated_proteins_fa
    mutated_proteins_annot = result.mutated_proteins_annot
}
