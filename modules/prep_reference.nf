
/*
 * PREP_REFERENCE
 *   - Build .fai and .dict for FASTA
 *   - Optionally derive "main chromosomes" FASTA (matching existing contig naming)
 *   - Build HISAT2 index (packed as tar.gz)
 *   - If GFF3 provided: filter it to the chosen contig naming, bgzip+tabix, make CDS BED
 *
 * Emits a single bundle channel with 7 elements:
 *   (ref_fa, ref_fa.fai, ref.dict, hisat2_index_tar, gff3_filtered_gz, gff3_filtered_gz.tbi, cds_bed)
 *
 * Notes:
 *   - If no GFF3 file provided, emits empty placeholders ('') for gff3 outputs and CDS BED.
 */

// Params 
def USE_MAIN_CHR    = params.use_main_chr    ?: false
def TOOLBOX_IMAGE   = params.toolbox_image   ?: 'ghcr.io/your-org/variant-caller-toolbox:latest'
def RGX_MAIN        = params.mainchr_regex   ?: '^(chr)?([0-9]{1,2}|X|Y|M|MT)$'  // tweak if needed

// Processes 

process FAIDX_DICT {
  tag "${ref_fa.simpleName}"
  label 'refprep'
  container TOOLBOX_IMAGE
  cpus 2
  memory '4 GB'
  publishDir "${params.outdir}/01_reference", mode: 'copy', overwrite: true

  input:
    path ref_fa

  output:
    tuple path("${base}.fasta"),
          path("${base}.fasta.fai"),
          path("${base}.dict")

  script:
    def base = ref_fa.baseName.replaceAll(/\.fa(sta)?$/, '')
    """
    set -euo pipefail
    cp -f ${ref_fa} ${base}.fasta

    if [ ! -s ${base}.fasta.fai ]; then
      samtools faidx ${base}.fasta
    fi

    if [ ! -s ${base}.dict ]; then
      if command -v gatk >/dev/null 2>&1; then
        gatk CreateSequenceDictionary -R ${base}.fasta -O ${base}.dict
      elif command -v picard >/dev/null 2>&1; then
        picard CreateSequenceDictionary R=${base}.fasta O=${base}.dict
      else
        echo "ERROR: Need GATK or Picard for CreateSequenceDictionary" >&2
        exit 2
      fi
    fi
    """
}

process MAKE_MAINCHR_FASTA {
  tag "${ref_fa.simpleName}.mainChr"
  label 'refprep'
  container TOOLBOX_IMAGE
  cpus 2
  memory '3 GB'
  publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

  input:
    tuple path(ref_fa), path(fai), path(dict)

  output:
    tuple path("${base}.fasta"),
          path("${base}.fasta.fai"),
          path("${base}.dict")

  when:
    USE_MAIN_CHR

  script:
    def base = ref_fa.baseName.replaceAll(/\.fa(sta)?$/, '') + ".mainChr"
    // Build a list of contigs matching the regex from the input FAI
    def awkCond = RGX_MAIN
    """
    set -euo pipefail

    awk 'BEGIN{OFS="\\t"} {print \$1}' ${fai} | \
      perl -ne 'print if /${awkCond}/' > keep.list

    # Extract selected contigs
    samtools faidx ${ref_fa} -r keep.list > ${base}.fasta

    samtools faidx ${base}.fasta

    if command -v gatk >/dev/null 2>&1; then
      gatk CreateSequenceDictionary -R ${base}.fasta -O ${base}.dict
    else
      picard CreateSequenceDictionary R=${base}.fasta O=${base}.dict
    fi
    """
}

process HISAT2_INDEX {
  tag "${ref_fa.simpleName}.hisat2"
  label 'refprep'
  container TOOLBOX_IMAGE
  cpus Math.max(4, (params.threads ?: 8) as int)
  memory '8 GB'
  publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

  input:
    path ref_fa

  output:
    path "${ref_fa.baseName}.hisat2.tar.gz"

  script:
    def base = ref_fa.baseName
    """
    set -euo pipefail

    if [ -s ${base}.hisat2.tar.gz ]; then
      echo "HISAT2 tarball exists. Skipping rebuild."
      exit 0
    fi

    hisat2-build -p ${task.cpus} ${ref_fa} ${base}.hisat2
    tar -czf ${base}.hisat2.tar.gz ${base}.hisat2.*.ht2
    """
}

process FILTER_GFF3_TO_REF {
  tag "${gff3.simpleName}->${fai.baseName}"
  label 'refprep'
  container TOOLBOX_IMAGE
  cpus 1
  memory '2 GB'
  publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

  input:
    path gff3
    path fai

  output:
    tuple path("${gff3.baseName}.filtered.gff3.gz"),
          path("${gff3.baseName}.filtered.gff3.gz.tbi")

  when:
    gff3
  script:
    """
    set -euo pipefail

    awk '{print \$1}' ${fai} > contigs.txt
    # Keep only features with seqid present in FASTA FAI, preserve directives/comments
    awk 'BEGIN{while((getline<"contigs.txt")>0) c[\$1]=1}
         /^#/ {print; next}
         c[\$1]==1 {print}' ${gff3} > filtered.gff3

    bgzip -f filtered.gff3
    tabix -p gff ${'filtered.gff3.gz'}
    mv filtered.gff3.gz ${gff3.baseName}.filtered.gff3.gz
    mv filtered.gff3.gz.tbi ${gff3.baseName}.filtered.gff3.gz.tbi
    """
}

process CDS_BED_FROM_GFF3 {
  tag "cds->bed"
  label 'refprep'
  container TOOLBOX_IMAGE
  cpus 1
  memory '1 GB'
  publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

  input:
    path gff3_gz

  output:
    path "${gff3_gz.baseName}.cds.bed"

  when:
    gff3_gz

  script:
    """
    set -euo pipefail
    # BED6 from CDS features. Start: 0-based.
    tabix -f -p gff ${gff3_gz} >/dev/null 2>&1 || true
    zcat ${gff3_gz} | awk -F'\\t' '
      BEGIN{OFS="\\t"}
      /^#/ {next}
      \$3=="CDS" {
        # parse attributes for ID or Parent or gene_name
        id="."
        n=split(\$9,a,";")
        for(i=1;i<=n;i++){
          split(a[i],b,"=")
          if(b[1]=="ID"){id=b[2]}
          else if(b[1]=="Parent" && id=="."){id=b[2]}
          else if(b[1]=="gene_name" && id=="."){id=b[2]}
        }
        start=\$4-1; if(start<0) start=0
        print \$1, start, \$5, id, 0, \$7
      }' > ${gff3_gz.baseName}.cds.bed
    """
}

// Wrapper workflow 

workflow prep_reference {

  take:
    ref_fasta          // required
    gff3_file          // optional (use Channel.fromPath(..., optional:true) in main)

  main:
    // 1) Base FASTA indices
    CH_BASE = FAIDX_DICT(ref_fasta)

    // 2) Choose reference: full or main-chr
    CH_CHOSEN = USE_MAIN_CHR ? MAKE_MAINCHR_FASTA(CH_BASE) : CH_BASE

    // 3) Build HISAT2 index tarball from chosen FASTA
    CH_HISAT = CH_CHOSEN.map { fa, fai, dict -> fa }.set { _fa_only }
    CH_HISAT = HISAT2_INDEX(_fa_only)

    // 4) If GFF3 provided: filter to chosen contigs, bgzip+tabix, make CDS BED
  gff3_file
    .filter { it }
    .combine( CH_CHOSEN.map { fa, fai, dict -> fai } )
    .filter { gff3, fai -> gff3 && gff3.toString().trim() }
    .set { GFF3_AND_FAI }

  def CH_GFF3_FILTERED = FILTER_GFF3_TO_REF(GFF3_AND_FAI)
  CH_GFF3 = CH_GFF3_FILTERED.ifEmpty { Channel.of( tuple('', '') ) }

  def CH_CDS_RAW = CDS_BED_FROM_GFF3(CH_GFF3_FILTERED.map{ g,t -> g })
  CH_CDS  = CH_CDS_RAW.ifEmpty { Channel.of( path('') ) }

    // 5) Emit tidy bundle
    REF_BUNDLE = CH_CHOSEN
      .combine( CH_HISAT )
      .combine( CH_GFF3 )
      .combine( CH_CDS )
      .map { fa, fai, dict,
             hisat_tar,
             gff3_filt, gff3_tbi,
             cds_bed ->
        tuple(fa, fai, dict, hisat_tar, gff3_filt, gff3_tbi, cds_bed)
      }

  emit:
    ref_bundle = REF_BUNDLE
}
