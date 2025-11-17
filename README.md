# Immunopeptidomics Variant Calling Workflow
Nextflow DSL2 pipeline for calling RNA-seq single-nucleotide variants (SNVs) and translating them into peptide sequences suitable for immunopeptidomics studies. The workflow trims raw reads, aligns them with HISAT2, performs variant calling (bcftools, GATK, or FreeBayes), annotates functional consequences, and extracts/ translates coding sequences carrying protein-altering variants.

## Key Capabilities
- Fastp-based read QC/trim, HISAT2 alignment, duplicate handling, and optional GATK SplitNCigarReads for RNA alignments.
- Variant calling with bcftools, GATK HaplotypeCaller, or FreeBayes, followed by bcftools consequence annotation.
- Automatic filtering for protein-altering variants, CDS extraction, and peptide translation for downstream immunopeptidomics.
- Flexible reference handling: build indexes/GFF3 derivatives on the fly or supply a pre-built bundle.
- Containerised (Docker/Apptainer) or Conda-based execution for reproducible environments.

## Repository Layout
- `main_SNV.nf` – entry-point pipeline defining orchestration and module wiring.
- `modules/` – reusable DSL2 modules (reference prep, alignment, callers, annotation, peptide generation).
- `nextflow.config` – default parameters, profiles, and resource settings.
- `env/base.yml` – Conda environment specification covering all required bioinformatics tools.
- `Dockerfile` – build recipe for a Micromamba-based toolbox image mirroring `env/base.yml`.
- `README.md` – project documentation (this file).

## Prerequisites
- Nextflow `>= 23.04` (DSL2 compatible).
- Either:
  - Container runtime (Docker, Podman, or Apptainer/Singularity), **or**
  - Conda/Mamba for creating environments from `env/base.yml`.
- Access to required input data: FASTQ reads, reference FASTA (+ optional GFF3/dbSNP).

## Getting Started
1. Clone the repository  
   ```bash
   git clone git@github.com:Lisa-G-88/variant_calling_immunopeptidogenomics.git
   cd variant_calling_immunopeptidogenomics
   ```
2. Choose an execution environment:
   - **Docker/Podman**  
     Build the toolbox image (if you do not already host it in a registry):
     ```bash
     docker build -t variant-caller-toolbox:latest .
     ```
     Launch Nextflow with `-profile docker` (default image name matches the config).
   - **Apptainer/Singularity**  
     Use `-profile docker` and ensure Apptainer is configured to run Docker images.
   - **Conda/Mamba**  
     ```bash
     mamba env create -f env/base.yml  # or conda env create
     conda activate variant-toolbox
     ```
     Launch Nextflow without the docker profile (or add a custom `-profile conda` if you extend the config).

## Running the Pipeline
Basic execution example:
```bash
nextflow run main_SNV.nf \
  --reads "data/*_{R1,R2}.fastq.gz" \
  --ref_fasta references/genome.fa \
  --gff3_file references/annotation.gff3.gz \
  --dbsnp_path references/dbsnp.vcf.gz \
  --outdir results \
  --caller bcftools \
  -profile docker
```

### Key Parameters
- `--reads` (required): FASTQ glob or pair pattern (`*_R{1,2}.fastq.gz`).
- `--ref_fasta` (required): reference genome FASTA.
- `--gff3_file`: GFF3 annotation (enables consequence annotation + CDS derivation).
- `--dbsnp_path`: dbSNP VCF or directory of VCFs for rsID annotation.
- `--caller`: `bcftools`, `gatk`, or `freebayes` (default: `bcftools`).
- `--outdir`: output directory root (`./results` by default).
- Reference bundle overrides (set `--prebuilt_bundle true` and provide `--ref_fai`, `--ref_dict`, `--hisat2_index`, `--gff3_filtered`, `--gff3_filtered_tbi`, `--cds_bed` to reuse precomputed assets).
- `--threads`, `--toolbox_image`, and other defaults are configured in `nextflow.config`.

A full list of parameters is documented inline in `main_SNV.nf` and `nextflow.config`.

## Outputs
Outputs are published beneath `--outdir`:
- `qc/fastp/` – HTML/JSON reports and trimmed FASTQs.
- `alignment/` – HISAT2 alignments and intermediate BAM/SAM files.
- `reference/` – Generated `.fasta`, `.fai`, `.dict`, HISAT2 indexes, filtered GFF3, and CDS BED (when reference prep runs).
- `vcf/` – Caller-specific VCFs, annotated variants, protein-altering subsets, CDS BEDs, and translated sequences.
- Additional logs and work directories follow Nextflow conventions.

## Reference Bundles
By default the workflow will build indexes and GFF3 derivatives. If you already maintain a reference bundle, supply the `--prebuilt_bundle true` flag along with the extra file paths so the pipeline skips rebuilding and reuses your assets.

## Development & Branching
- `main` – primary pipeline code and modules (recommended for production).
- `support` – optional branch for staging enhancements or exploratory work.

Feel free to extend the workflow with additional modules (e.g., alternative aligners) or new profiles. Contributions and issues are welcome.
