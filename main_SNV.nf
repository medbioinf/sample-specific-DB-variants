#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Params

params.reads               = null           // e.g. "data/*_{R1,R2}.fastq.gz" or "data/*.fastq.gz"
params.ref_fasta           = null           // path to FASTA
params.dbsnp_path          = null           // e.g., ".vcf.gz" or dir with multiple .vcf.gz files
params.chrom_rename_map    = null           // optional 2-col map for --rename-chrs
params.outdir              = "./results"

// Optional: provide GFF3 so prep_reference can filter + make CDS BED
params.gff3_file           = null
params.use_main_chr        = true

// Prebuilt bundle (set true only if you supply all below)
params.prebuilt_bundle     = false
params.ref_fai             = null           // .fai
params.ref_dict            = null           // .dict
params.hisat2_index        = null           // dir with *.ht2 OR a .tar.gz
params.gff3_filtered       = null           
params.gff3_filtered_tbi   = null           
params.cds_bed             = null            

// General
params.threads             = 8
params.toolbox_image       = 'ghcr.io/your-org/variant-caller-toolbox:latest'

// Caller: bcftools | gatk | freebayes
params.caller              = 'bcftools'    // bcftools | gatk | freebayes
params.min_dp              = 3             // minimum depth of coverage (total number of reads (REF+ALT))
params.min_alt             = 3             // minimum ALT-supporting reads
params.protein_mode        = 'merged'      // merged | per_variant


// Modules
include { prep_reference } from './modules/prep_reference.nf'
include { prep_reads     } from './modules/prep_reads.nf'
include { align_hisat2   } from './modules/align_hisat2.nf'  
include { sort_index_bam } from './modules/sort_index_bam.nf'
include { fixmate        } from './modules/fixmate.nf'
include { markdup        } from './modules/markdup.nf'
include { split_n_cigar  } from './modules/split_n_cigar.nf'
include { gatk_haplotypecaller } from './modules/callers/gatk_haplotypecaller.nf'
include { bcftools_call  } from './modules/callers/bcftools_call.nf'
include { freebayes_call } from './modules/callers/freebayes_call.nf'
include { annotate_variants } from './modules/annotate_variants.nf'
include { filter_and_extract } from './modules/filter_and_extract.nf'
include { extract_cds_bed } from './modules/extract_cds_bed.nf'
include { extract_sequences } from './modules/extract_sequences.nf'
include { translate_proteins } from './modules/translate_proteins.nf'
include { per_variant_proteins } from './modules/per_variant_proteins.nf'
include { match_rna_editing } from './modules/match_rna_editing.nf'
include { final_protein_db } from './modules/final_protein_db.nf'


// Validation

if (!params.reads)     exit 1, "ERROR: --reads is required"
if (!params.ref_fasta) exit 1, "ERROR: --ref_fasta is required"
if (!params.outdir)    exit 1, "ERROR: --outdir is required"
if (params.caller?.toLowerCase() !in ['bcftools', 'gatk', 'freebayes']) {
  exit 1, "ERROR: --caller must be one of: bcftools | gatk | freebayes"
}

// Helper

def havePrebuilt() {
  (params.prebuilt_bundle as boolean) &&
  params.ref_fai && params.ref_dict && params.hisat2_index &&
  (
    (!params.gff3_filtered && !params.gff3_filtered_tbi) ||
    (params.gff3_filtered && params.gff3_filtered_tbi)
  )
}

// --------------
// workflow
// --------------
workflow {

  def caller  = params.caller?.toLowerCase()
  def useGatk = (caller == 'gatk')

  // inputs
  // simple detector for paired pattern
  def isPairedPattern = (params.reads.contains('{R1,R2}') || params.reads.contains('{1,2}'))

  ch_reads = isPairedPattern
    ? Channel.fromFilePairs(params.reads, checkIfExists: true)
             .map { id, files -> tuple(id as String, files) }
    : Channel.fromPath(params.reads, checkIfExists: true)
             .map { p -> tuple(p.baseName.replaceAll(/\.fastq|fq(\.gz)?$/, ''), p) }

  ch_ref_fasta = Channel.fromPath(params.ref_fasta, checkIfExists: true)
  ch_gff3      = params.gff3_file
    ? Channel.fromPath(params.gff3_file, checkIfExists: true)
    : Channel.empty()

  ch_dbsnp  = params.dbsnp_path       ? Channel.fromPath(params.dbsnp_path, checkIfExists: true)
                                    : Channel.value(null)
  ch_chrmap = params.chrom_rename_map ? Channel.fromPath(params.chrom_rename_map, checkIfExists: true)
                                    : Channel.value(null)

  // build or use prebuilt bundle
  def ch_ref_bundle   // declare a channel variable we'll assign

  if (havePrebuilt()) {
    def ref_fa   = file(params.ref_fasta)
    def ref_fai  = file(params.ref_fai)
    def ref_dict = file(params.ref_dict)
    def hisat2   = file(params.hisat2_index)

    // only wrap with file() when they exist; otherwise use null
    def gff3_f   = params.gff3_filtered     ? file(params.gff3_filtered)     : null
    def gff3_tbi = params.gff3_filtered_tbi ? file(params.gff3_filtered_tbi) : null
    def cds_bed  = params.cds_bed           ? file(params.cds_bed)           : null

    // create a single-emission channel carrying the 7-tuple
    ch_ref_bundle = Channel.value( tuple(ref_fa, ref_fai, ref_dict, hisat2, gff3_f, gff3_tbi, cds_bed) )

    log.info "Reference bundle is used"
  }
  else {
    // prep_reference must define: out: tuple val(fa), val(fai), val(dict), path(hisat2_tar_or_dir), path(gff3_f), path(gff3_tbi), path(cds_bed)
    prep_reference(ch_ref_fasta, ch_gff3)
    ch_ref_bundle = prep_reference.out.ref_bundle
      .map { fa, fai, dict, hisat2, gff, gff_tbi, cds ->
        log.info "Refercne bundle is created"
        tuple(fa, fai, dict, hisat2, gff, gff_tbi, cds)
      }

    log.info "Building reference bundle via prep_reference."
  }

  // split for downstream modules
  ch_ref_for_align = ch_ref_bundle.map { fa, fai, dict, hisat2, gff, gff_tbi, cds -> hisat2 }
  ch_ref_for_call  = ch_ref_bundle.map { fa, fai, dict, hisat2, gff, gff_tbi, cds -> tuple(fa, dict) }
  ch_annot_res     = ch_ref_bundle.map { fa, fai, dict, hisat2, gff, gff_tbi, cds -> tuple(gff, gff_tbi, cds) }
  
  // onl fasta is used for bcftools and freebaes
  def ref_fa_only  = ch_ref_bundle.map { fa, fai, dict, hisat2, gff, gff_tbi, cds -> fa } // fasta-only channel 

  //-------------------------------
  // Variant Calling workflow steps
  //-------------------------------

  // 1) fastq QC/trim with fastp
  prep_reads(ch_reads)    // emits (sample_id, *.trim.fastq.gz)
  def ch_trimmed = prep_reads.out.trimmed

  // 2) align with HISAT2 (module accepts index dir or .tar.gz)
  align_hisat2( ch_trimmed, ch_ref_for_align )  // emits (sample_id, ${sample_id}.sam) or (sample_id, ${sample_id}.unsorted.bam)
  def ch_align_out = align_hisat2.out[0]
  ch_align_out.subscribe { tuple ->
    def sid = tuple instanceof java.util.List ? tuple[0] : tuple?.getAt(0)
    if (sid) log.info "Alignment is complited (${sid})"
  }
  def ch_sam_or_bam = ch_align_out

  // 3) sort + index BAM
  sort_index_bam(ch_sam_or_bam)    // emits (sample_id, ${sample}.sorted.bam) + .bai
  def ch_sorted_bam = sort_index_bam.out.bam

  // attachment of mate tages (for PE-reads only)
  def ch_for_markdup  = isPairedPattern ? fixmate(ch_sorted_bam) : ch_sorted_bam
  
  // 4) mark duplicates
  markdup(ch_for_markdup)   // emits (sample_id, ${sample}.mark
  def ch_dedup_bam = markdup.out.bam

  // 5) split N CIGAR reads if GATK will be used 
  def ch_ref_for_gatk = useGatk ? ch_ref_bundle.map { fa, fai, dict, hisat2, gff, gff_tbi, cds -> tuple(fa, fai, dict) } : null
  def ch_split = useGatk ? split_n_cigar(ch_dedup_bam, ch_ref_for_gatk) : null   // emits (sample_id, *.split.bam) when used
  
  def bam_for_calling = useGatk
    ? ch_split.split_bam      // emits (sample_id, *.split.bam)
    : ch_dedup_bam            // emits (sample_id, *.markdup.bam)

  // 6) call variants: bcftools | gatk | freebayes 
  def raw_vcf =
  (caller == 'bcftools') ? bcftools_call(bam_for_calling, ref_fa_only, params.min_dp as int, params.min_alt as int) :
  (useGatk)              ? gatk_haplotypecaller(bam_for_calling, ch_ref_for_gatk, params.min_dp as int, params.min_alt as int) :
  (caller == 'freebayes')? freebayes_call(bam_for_calling, ref_fa_only, params.min_dp as int, params.min_alt as int) :
  { exit 1, "Unknown --caller '${params.caller}'. Use: bcftools | gatk | freebayes" }()

  // take the unified outlet for downstream
  def ch_vcf_channel = raw_vcf.vcf_gz
  ch_vcf_channel.subscribe { sid, vcf ->
    if (sid) log.info "Variants calling step is completed (${sid})"
  }
  def ch_vcf_gz = ch_vcf_channel
  // combine filtered .vcf + index into one channel
  def vcf_tbi_ch = ch_vcf_gz.map { sid, vcf -> tuple(sid, vcf, file("${vcf}.tbi"))}  // emits (sample_id, *.vcf.gz, *.vcf.gz.tbi)
 
  //--------------------------------------
  // Peptide Identification workflow steps
  //--------------------------------------

  // Reference FASTA (.fa/.fai)
  def ref_pair_for_ann = ch_ref_bundle
    .map { fa, fai, dict, hisat2, gff, gff_tbi, cds -> tuple(fa, fai) }
    .filter { fa, fai -> fa && fai }
    .ifEmpty { exit 1, "Reference FASTA (.fa) or index (.fai) missing." }

  def ref_pair_for_cds = ch_ref_bundle
    .map { fa, fai, dict, hisat2, gff, gff_tbi, cds -> tuple(fa, fai) }
    .filter { fa, fai -> fa && fai }
    .ifEmpty { exit 1, "Reference FASTA (.fa) or index (.fai) missing." }

  // GFF3 (.gff3.gz/.tbi) (and optional CDS BED)
  def gff3_for_ann = ch_annot_res
    .map { gff, gff_tbi, cds -> tuple(gff, gff_tbi, cds) }
    .filter { gff, gff_tbi, cds -> gff && gff_tbi }
    .ifEmpty { exit 1, "GFF3 (.gff3.gz) or its index (.tbi) missing. Provide --gff3_filtered/--gff3_filtered_tbi (or supply --gff3_file so prep_reference can build them)." }

  def gff3_for_seq = ch_annot_res
    .map { gff, gff_tbi, cds -> tuple(gff, gff_tbi, cds) }
    .filter { gff, gff_tbi, cds -> gff && gff_tbi }
    .ifEmpty { exit 1, "GFF3 (.gff3.gz) or its index (.tbi) missing. Provide --gff3_filtered/--gff3_filtered_tbi (or supply --gff3_file so prep_reference can build them)." }

  // combine inputs for annotation: filtered .vcf/.tbi + ref.fa/.fai + gff3.gz/.tbi
  def ann_input_ch = vcf_tbi_ch
    .combine(ref_pair_for_ann)
    .combine(gff3_for_ann)
    .map { sample_id, filtered_vcf, filtered_vcf_tbi, ref_fa, ref_fai, gff3_path, gff3_tbi_path, cds_bed_ignore ->
      tuple(sample_id, filtered_vcf, filtered_vcf_tbi, ref_fa, ref_fai, gff3_path, gff3_tbi_path)
    }
  // 1) annotate variants
  annotated_vcf = annotate_variants(ann_input_ch)   // emits annotated_vcf

  // 2) filter for protein-altering variants and extracts affeected trasncripts
  def protein_altering_channel = filter_and_extract(annotated_vcf)    // emits protein_altering_bundle: tuple(sid, prot_alter.vcf.gz, .tbi, positions.tsv, bed, affected_transcripts.txt)
  protein_altering_channel.subscribe { tuple ->
    def sid = tuple instanceof java.util.List ? tuple[0] : tuple?.getAt(0)
    if (sid) log.info "Protein-altering variants are identified (${sid})"
  }
  def protein_altering_bundle = protein_altering_channel
  
  def filtered_vcf_all = protein_altering_bundle.map { sid, vcf, tbi, pos, bed, tx -> tuple(sid, vcf, tbi) }
  def transcripts      = protein_altering_bundle.map { sid, vcf, tbi, pos, bed, tx -> tuple(sid, pos, bed, vcf, tbi, tx) }
  def proteinMode      = (params.protein_mode ?: 'merged')?.toLowerCase()
  def useMergedMode    = (proteinMode == 'merged')
  
  // Optional step: match RNA editing sites (A to I RNA editing events)
  def editing_catalog   = params.rna_editing ? file(params.rna_editing) : null
  def editing_matches_ch = null
  if (editing_catalog) {
    def editing_input = protein_altering_bundle.map { sid, vcf, tbi, pos, bed, tx ->
      tuple(sid, vcf, tbi, editing_catalog)
    }
    match_rna_editing(editing_input)
    editing_matches_ch = match_rna_editing.out.editing_matches
  }
  
  // 3) extract BED of CDS regions for filtered transcripts
  sequences_input_ch = transcripts
    .combine(gff3_for_seq)
    .map { sample_id, positions_tsv, prot_alt_bed, filtered_vcf, filtered_vcf_tbi, affected_tx, gff3_path, gff3_tbi_path, cds_bed_ignore ->
      tuple(sample_id, positions_tsv, prot_alt_bed, filtered_vcf, filtered_vcf_tbi, affected_tx, file(gff3_path), file(gff3_tbi_path))
    }
  cds_bed = extract_cds_bed(sequences_input_ch)  // emits cds_affected_bed: tuple(sid, positions.tsv, prot_alt.bed, prot_alt.vcf.gz, .tbi, affected_tx.txt, cds_affected.bed)

  // 4/5) translate CDS sequences to proteins based on mode
  def mutated_proteins_fa
  def mutated_proteins_annot

  if (useMergedMode) {
    // extract FASTA sequences for those CDS regions (all variants together)
    def cds_input_ch = cds_bed
      .combine(ref_pair_for_cds)
      .map { sid, positions_tsv, prot_alt_bed, filtered_vcf, filtered_vcf_tbi, affected_tx, cds_affected_bed, ref_fa, ref_fai ->
        tuple(sid, file(ref_fa), file(ref_fai), cds_affected_bed, filtered_vcf, filtered_vcf_tbi)
      }
    cds_sequences = extract_sequences(cds_input_ch)

    // translate merged CDS per transcript (existing behavior)
    def vcf_for_translate = protein_altering_bundle.map { sid, vcf, tbi, pos, bed, tx -> tuple(sid, vcf, tbi) }
    def mutated_cds_ch = cds_sequences.mutated_cds_fa
      .map { f -> tuple(f.baseName.replace('.mutated_cds',''), f) }   // (sid, mutated_cds.fa)
    def original_cds_ch = cds_sequences.original_cds_fa
      .map { f -> tuple(f.baseName.replace('.original_cds',''), f) }  // (sid, original_cds.fa)
    def translate_input_ch_base = mutated_cds_ch
      .join(original_cds_ch)                                          // [sid, mutated_cds.fa, original_cds.fa]
      .join(vcf_for_translate)                                        // [sid, mutated_cds.fa, original_cds.fa, vcf, tbi]
    def translate_input_ch
    if (editing_catalog) {
      translate_input_ch = translate_input_ch_base
        .join(editing_matches_ch)
        .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5]) }
    }
    else {
      def placeholder_path = "${workflow.projectDir}/.empty_rna_editing.txt"
      def placeholder_file = new File(placeholder_path)
      if (!placeholder_file.exists()) {
        placeholder_file.text = ''
      }
      def placeholder_nf_file = file(placeholder_path)
      translate_input_ch = translate_input_ch_base
        .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], placeholder_nf_file) }
    }
    translate_proteins(translate_input_ch)

    mutated_proteins_fa    = translate_proteins.out.mutated_proteins_fa
      .map { path_obj ->
        def base = path_obj?.baseName ?: path_obj?.name ?: path_obj?.toString()
        def sid = base?.replace('.mutated_proteins','')
        log.info "CDSs are translated to proteins (${sid ?: 'unknown sample'})"
        path_obj
      }
    mutated_proteins_annot = translate_proteins.out.mutated_proteins_annot
  }
  else {
    // per-variant mode: apply variants one-by-one and translate
    def per_variant_input = cds_bed
      .combine(ref_pair_for_cds)
      .map { sid, positions_tsv, prot_alt_bed, filtered_vcf, filtered_vcf_tbi, affected_tx, cds_affected_bed, ref_fa, ref_fai ->
        tuple(sid, file(ref_fa), file(ref_fai), filtered_vcf, filtered_vcf_tbi, cds_affected_bed)
      }
    per_variant_proteins(per_variant_input)
    mutated_proteins_fa    = per_variant_proteins.out.mutated_proteins_fa
    mutated_proteins_annot = per_variant_proteins.out.mutated_proteins_annot
  }

  // 6) concatenate mutated proteins with: canonical database + decoy database (if provided)

  def canonical_db = params.canonical_proteins ? file(params.canonical_proteins) : null
  def decoy_db     = params.decoy_proteins ? file(params.decoy_proteins) : null

  if (canonical_db && decoy_db) {
    def mutated_channel = mutated_proteins_fa.map { fa ->
      def base = fa?.baseName ?: fa?.name ?: fa?.toString()
      def sid = base?.replace('.mutated_proteins','')?.replace('.per_variant_proteins','')
      tuple(sid, fa)
    }
    def annot_channel   = mutated_proteins_annot.map { annot ->
      def base = annot?.baseName ?: annot?.name ?: annot?.toString()
      def sid = base?.replace('.mutated_proteins.annot','')?.replace('.per_variant_proteins.annot','')
      tuple(sid, annot)
    }
    def final_db_input  = mutated_channel
      .join(annot_channel)
      .map { row -> tuple(row[0], row[1], row[2], canonical_db, decoy_db) }
    final_protein_db(final_db_input)
    final_protein_db.out.extended_protein_db.subscribe { tuple ->
      def sid = tuple instanceof java.util.List ? tuple[0] : tuple?.getAt(0)
      if (sid) log.info "Canincal protein database is extended with mutated proteins (${sid})"
    }
  }
  else if (canonical_db || decoy_db) {
    log.warn "Both --canonical_proteins and --decoy_proteins are required to build the extended FASTA; skipping final concatenation."
  }
  else {
    log.info "Canonical/decoy FASTA inputs not provided; skipping extended protein database concatenation."
  }
  
}
