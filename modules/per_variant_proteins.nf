nextflow.enable.dsl=2
// Per-variant protein generation: one protein per variant x transcript.
process per_variant_proteins {
    tag "$sample_id"
    publishDir "${params.outdir}/07_translate_proteins", mode: 'copy', overwrite: true, dynamic: true
    cpus { params.threads ?: 8 }
    container (params.toolbox_image ?: 'ghcr.io/your-org/variant-caller-toolbox:latest')

    input:
    tuple val(sample_id),
         path(reference_fasta),
         path(reference_fai),
         path(protein_vcf),
         path(protein_vcf_tbi),
         path(cds_bed)

    when:
    reference_fasta.exists() && reference_fai.exists() &&
    protein_vcf.exists() && protein_vcf_tbi.exists() &&
    cds_bed.exists()

    output:
    path("${sample_id}.per_variant_proteins.fa"),        emit: mutated_proteins_fa
    path("${sample_id}.per_variant_proteins.annot.txt"), emit: mutated_proteins_annot

    script:
    """
    set -euo pipefail

    # No variants? emit empty files
    if [ ! -s "${protein_vcf}" ]; then
      : > ${sample_id}.per_variant_proteins.fa
      : > ${sample_id}.per_variant_proteins.annot.txt
      exit 0
    fi

    # 1) atomize/normalize
    bcftools norm --atomize -m -both -f ${reference_fasta} ${protein_vcf} -Oz -o norm.vcf.gz
    bcftools index -f -t norm.vcf.gz

    export SAMPLE_ID="${sample_id}"
    export CALLER="${params.caller}"

    # 2) build variants.tsv (id, chrom, pos, ref, alt, tx list)
    python3 - <<'PY'
import os, re, subprocess, csv
bcsq_re = re.compile(r'BCSQ=([^;]+)')
cmd = ["bcftools", "view", "-H", "norm.vcf.gz"]
with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc, open("variants.tsv", "w", newline="") as out_handle:
    writer = csv.writer(out_handle, delimiter="\\t")
    idx = 0
    for raw in proc.stdout:
        raw = raw.rstrip('\\n')
        if not raw:
            continue
        parts = raw.split('\\t')
        if len(parts) < 8:
            continue
        chrom, pos, _id, ref, alt, _qual, _filt, info = parts[:8]
        match = bcsq_re.search(info)
        if not match:
            continue
        tx_ids = set()
        for entry in match.group(1).split(','):
            bits = entry.split('|')
            for b in bits:
                if b.startswith('ENST'):
                    tx_ids.add(b)
                elif b.startswith('transcript:'):
                    t = b.split(':', 1)[1]
                    if t.startswith('ENST'):
                        tx_ids.add(t)
        if not tx_ids:
            continue
        idx += 1
        writer.writerow([idx, chrom, pos, ref, alt, ','.join(sorted(tx_ids))])
if not os.path.exists('variants.tsv') or os.path.getsize('variants.tsv') == 0:
    open('variants.tsv', 'w').close()
PY

    if [ ! -s variants.tsv ]; then
      echo "No variants with transcript annotations found."
      : > ${sample_id}.per_variant_proteins.fa
      : > ${sample_id}.per_variant_proteins.annot.txt
      exit 0
    fi

    : > ${sample_id}.per_variant.merged_cds.fa
    : > ${sample_id}.per_variant_proteins.annot.txt
    : > ${sample_id}.per_variant_proteins.fa

    python3 - <<'PY'
import csv, os, shlex, subprocess
from pathlib import Path

variant_file = 'variants.tsv'
cds_bed = '${cds_bed}'
ref_fasta = '${reference_fasta}'
annot_out = '${sample_id}.per_variant_proteins.annot.txt'
cds_out = '${sample_id}.per_variant.merged_cds.fa'
sample = os.environ.get("SAMPLE_ID", "sample")
caller = os.environ.get("CALLER", "")

tx_strand_map = {}


def run(cmd):
    subprocess.run(cmd, check=True)

def get_aa_labels_and_genes(vcf_path, tx):
    aa_labels = []
    genes = []
    utr_flag = False
    try:
        res = subprocess.run(
            ['bcftools', 'query', '-f', '%INFO/BCSQ\\n', vcf_path],
            capture_output=True, text=True, check=True
        )
        for line in res.stdout.splitlines():
            for entry in line.split(','):
                bits = entry.split('|')
                tx_match = (tx in bits) or (f'transcript:{tx}' in bits)
                if not tx_match:
                    continue
                gene = bits[1] if len(bits) > 1 else ''
                if gene and gene not in genes:
                    genes.append(gene)
                consequence = bits[0] if bits else ''
                aa = bits[5] if len(bits) > 5 else ''
                if aa and aa not in {'.', '-'} and aa not in aa_labels:
                    aa_labels.append(aa)
                if 'utr' in consequence.lower():
                    utr_flag = True
    except Exception:
        pass
    return aa_labels, genes, utr_flag

def translate(seq, table_id):
    seq = seq.upper().replace('U','T')
    table1 = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    }
    table2 = dict(table1)
    table2.update({'AGA':'*','AGG':'*','ATA':'M','TGA':'W'})
    table = table2 if table_id == 2 else table1
    prot = []
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        prot.append(table.get(codon, 'X'))
    return ''.join(prot)

bed_cache = Path("bed_cache")
bed_cache.mkdir(exist_ok=True)

# cache CDS intervals for each transcript so we don't rewrite bed files
cds_map = {}
with open(cds_bed) as bed_in:
    for line in bed_in:
        if not line.strip():
            continue
        parts = line.rstrip('\\n').split('\\t')
        if len(parts) < 4:
            continue
        chrom, start, end, tx = parts[:4]
        strand = parts[5] if len(parts) > 5 else '+'
        phase = 0
        if len(parts) > 6:
            try:
                phase = int(parts[6])
            except Exception:
                phase = 0
        cds_map.setdefault(tx, []).append((chrom, int(start), int(end), strand, tx, phase))
        tx_strand_map[tx] = strand

def write_tx_bed(tx):
    bed_path = bed_cache / f"{tx}.bed"
    # always rebuild to incorporate any updated phase/ordering logic
    if bed_path.exists():
        bed_path.unlink()
    entries = cds_map.get(tx)
    if not entries:
        return None
    strand = entries[0][3]
    # For minus strand, maintain descending genomic order so concatenation matches transcript 5'->3'
    entries_sorted = sorted(entries, key=lambda x: x[1], reverse=(strand == '-'))
    with open(bed_path, 'w') as out:
        for chrom, start, end, strand, name, phase in entries_sorted:
            # Apply phase to stay in-frame: plus strand -> start+phase; minus strand -> end-phase
            if strand == '-':
                adj_start = start
                adj_end = end - phase
            else:
                adj_start = start + phase
                adj_end = end
            if adj_end <= adj_start:
                continue
            out.write(f"{chrom}\\t{adj_start}\\t{adj_end}\\t{name}\\t0\\t{strand}\\n")
    return bed_path

def revcomp(seq):
    comp = str.maketrans('ACGTNacgtn','TGCANtgcan')
    return seq.translate(comp)[::-1]

def merge_fasta(path, strand='+'):
    # Concatenate sequences from a multi-FASTA.
    # For minus strand, reverse-complement each record individually while preserving record order.
    seqs = []
    with open(path) as handle:
        cur = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur:
                    seqs.append(''.join(cur))
                    cur = []
            else:
                cur.append(line)
        if cur:
            seqs.append(''.join(cur))
    if strand == '-':
        seqs = [revcomp(s) for s in seqs]
    return ''.join(seqs)

seq_counter = {}
ref_seq_cache = {}

with open(variant_file) as vf:
    reader = csv.reader(vf, delimiter='\\t')
    for row in reader:
        if len(row) < 6:
            continue
        vid, chrom, pos, ref, alt, txs_str = row
        if not txs_str:
            continue

        # build per-variant VCF then bgzip/tabix for consensus
        with open('var.vcf', 'w') as vout:
            subprocess.run(['bcftools', 'view', '-h', 'norm.vcf.gz'], check=True, stdout=vout)
            subprocess.run(['bcftools', 'view', '-H', '-r', f'{chrom}:{pos}-{pos}', 'norm.vcf.gz'], check=True, stdout=vout)
        run(['bgzip', '-f', 'var.vcf'])
        run(['tabix', '-f', 'var.vcf.gz'])

        txs = [t for t in txs_str.split(',') if t]
        gc = ${params.genetic_code ?: 1}
        if chrom in ['MT','chrM','M']:
            gc = 2
        for tx in txs:
            seq_counter[tx] = seq_counter.get(tx, 0) + 1
            suffix = f'_{seq_counter[tx]}'

            bed_path = write_tx_bed(tx)
            if bed_path is None:
                continue

            alt_fa = f"{tx}.alt.fa"
            # Convert BED (0-based) to samtools region strings (1-based) so consensus can map variants.
            # Apply phase correctly per strand to avoid frame shifts.
            region_list = f"{tx}.regions.txt"
            with open(bed_path) as bed_in, open(region_list, 'w') as reg_out:
                for line in bed_in:
                    if not line.strip():
                        continue
                    cols = line.rstrip('\\n').split('\\t')
                    if len(cols) < 6:
                        continue
                    chrom, start, end, _name, _score, strand = cols[:6]
                    try:
                        start_i = int(start)
                        end_i = int(end)
                    except ValueError:
                        continue
                    if end_i <= start_i:
                        continue
                    reg_out.write(f"{chrom}:{start_i+1}-{end_i}\\n")

            run([
                'bash', '-c',
                f"samtools faidx {shlex.quote(ref_fasta)} -r {shlex.quote(region_list)} | "
                f"bcftools consensus var.vcf.gz -o {shlex.quote(alt_fa)}"
            ])
            tx_strand = tx_strand_map.get(tx, '+')
            merged_seq = merge_fasta(alt_fa, strand=tx_strand)
            if not merged_seq:
                continue
            merged_seq = ''.join(base if base.upper() in {'A','C','G','T','N'} else 'N' for base in merged_seq)
            # trim to codon boundary to avoid downstream frame drift
            rem = len(merged_seq) % 3
            if rem:
                merged_seq = merged_seq[:-rem]

            variant_label = f'{tx}_var{suffix}'
            variant_locus = f'{chrom}:{pos} {ref}>{alt}'
            aa_labels, genes, utr_flag = get_aa_labels_and_genes('var.vcf.gz', tx)

            if tx not in ref_seq_cache:
                ref_seq_cache[tx] = ''
                res = subprocess.run(
                    ['bedtools', 'getfasta', '-fi', ref_fasta, '-bed', str(bed_path), '-s', '-name'],
                    capture_output=True, text=True
                )
                if res.returncode == 0:
                    ref_seq_cache[tx] = ''.join(
                        line.strip()
                        for line in res.stdout.splitlines()
                        if line and not line.startswith('>')
                    )
            ref_seq = ref_seq_cache.get(tx, '')
            prot_ref = translate(ref_seq, gc) if ref_seq else ''
            prot_alt = translate(merged_seq, gc) if merged_seq else ''

            diffs = []
            clean_ref = prot_ref.replace('*', '')
            clean_alt = prot_alt.replace('*', '')
            max_len = max(len(clean_ref), len(clean_alt))
            for i in range(max_len):
                a = clean_ref[i] if i < len(clean_ref) else '-'
                b = clean_alt[i] if i < len(clean_alt) else '-'
                if a != b:
                    diffs.append(f"{i+1}{a}>{b}")
            if len(clean_ref) != len(clean_alt):
                diffs.append(f"len{len(clean_ref)}->{len(clean_alt)}")
            if diffs:
                aa_diff_str = ';'.join(diffs)
            else:
                aa_diff_str = 'variant hits the UTR - no amino acid change' if utr_flag else 'no amino acid change (CDS unchanged)'

            prot_path = f\"{sample}.per_variant_proteins.fa\"
            with open(cds_out, 'a') as out_handle, open(annot_out, 'a') as ann_handle, open(prot_path, 'a') as prot_handle:
                out_handle.write(f'>{tx}{suffix}\\n{merged_seq}\\n')
                ann_handle.write(f'>{tx}{suffix}\\n')
                ann_handle.write(f'variant_index: {suffix.lstrip(\"_\")}\\n')
                ann_handle.write(f'variant: {variant_label}\\n')
                ann_handle.write(f'variant_locus: {variant_locus}\\n')
                gene_str = ';'.join(genes) if genes else 'NA'
                ann_handle.write(f'amino_acid_change: {aa_diff_str}\\n')
                ann_handle.write(f'gene: {gene_str}\\n')
                ann_handle.write(f'caller: {caller}\\n')
                prot_handle.write(f'>{tx}{suffix}\\n{clean_alt}\\n')

print('Finished per-variant build.')
PY

    if [ ! -s ${sample_id}.per_variant.merged_cds.fa ]; then
      echo "No CDS sequences generated."
      : > ${sample_id}.per_variant_proteins.fa
      exit 0
    fi
    """
}

workflow {
  take:
    input_ch

  main:
    res = per_variant_proteins(input_ch)

  emit:
    mutated_proteins_fa    = res.mutated_proteins_fa
    mutated_proteins_annot = res.mutated_proteins_annot
}
