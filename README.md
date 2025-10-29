# Sample-Specific Variant Calling and Peptide Database Generation Workflow for Immunopeptidogenomics 

[Nextflow](https://www.nextflow.io/) DSL2 pipeline for calling RNA-seq variants and translating them into peptide sequences suitable for immunopeptidomics studies. The workflow trims raw reads, aligns them with HISAT2, performs variant calling (bcftools, GATK, or FreeBayes), annotates functional consequences, and extracts/ translates coding sequences carrying protein-altering variants.


<details>
  <summary><b>Introduction</b></summary>

### Biological background

Immunopeptidogenomics is a multidisciplinary field that combines genomics, transcriptomics, and immunopeptidomics to investigate how genetic variations influence the range of peptides displayed on the cell surface by [major histocompatibility complex (MHC)](https://www.sciencedirect.com/science/article/pii/S0091674909028371?via%3Dihub) proteins. Implementation of mass spectrometry–based peptide identification along with next-generation sequencing (NGS) data allows for linking distinct genetic variants to their potential immunological impacts.

Conventional proteomic studies rely on reference protein databases (e.g., [UniProt](https://www.uniprot.org/) or [Ensembl](https://www.ensembl.org/index.html)) representing the canonical gene products of human DNA. However, these standard databases do not account for sample-specific genetic modifications, such as single nucleotide variants (SNVs), deletions/inserts and RNA-editing events, which could lead to altered or non-canonical peptides. These variant peptides often go unidentified in standard analysis procedures, even though they are biologically relevant.

---

### Variant-derived peptides

Detection and characterization of variant-derived peptides is essential for the understanding of how mutations impact immune recognition. In cancer and infectious disease, such peptides can serve as neoantigens—novel, non-self antigens presented by either tumor or infected cells that are recognized by the immune system as foreign. Identifying these disease-associated peptides is informative for the development of immunotherapy and precise medical treatment.

From a systems biology perspective, building sample-specific protein databases with detected variants provides a more thorough and representative personalized view of antigen presentation. This allows for greater proteome coverage and improves upon the sensitivity of peptide identification in mass spectrometry workflows as well as highlighting previously obscured regions of the "dark immunopeptidome"— antigenic peptides not represented by standard reference databases.

---

### Bioinformatics background

The majority of workflows already established for immunopeptidogenomics workflows are tied to tools or specific environments. This constraining factor hinders flexibility, scalability, and overall workflow reproducibility. This workflow proposes an alternative approach by integrating a modular, containerized, and reproducible workflow based on Nextflow DSL2. The pipeline supports automated variant calling from RNA-seq data with user-defined parameters and integrates transcriptomic variants consisitng of nonsynonymous SNVs and indels with canonical proteomes through the creation of sample-specific protein databases. These protein databases allow for greater accuracy and transparency when identifying HLA-presented peptides downstream, and allow for easy adaptation of more workflows, even when computational requirements differ.

</details>


<details>
  <summary><b>Purpose & Features</b></summary>
<pre>
Enables automated discovery of genetic variants and creation of variant-containing peptide database.

### Main features:
  - Modular Nextflow DSL2 structure for transparency and flexibility
  - Support for three variant callers (user-defined choice)
  - Built-in Docker and Conda compatibility for reproducibility
</pre>
</details>

<details>
  <summary><b>Pipeline Overview</b></summary>
<pre>
1. Input Preparation
    -  Normalization of RNA-seq data; 
    -  Pre-processing of refernce and annotation data
2. Read Alignment
    -  Aligning reads to the reference genome
3. Post-processing
    -  Refining alignment results
4. Variant-calling
    -  Identifying variants using:BCFtools, FreeBayes or GATK
5. Variants Annotation
6. Extraction of Coding Sequences 
7. Generation of Variant-containing Sequences
8. Translation and Database Creation
    -  Converting nucleotide sequences into amino acid sequences 

</pre>
</details>

<details>
  <summary><b>Repository structure</b></summary>
<pre>
immunopeptidogenomics-workflow/
├── main_SNV.nf                   # Main Nextflow pipeline script
├── modules/                      # Individual DSL2 modules
│   ├── prep_reference.nf
│   ├── alignment_hisat2.nf
│   ├── markdup.nf
│   ├── split_n_cigar.nf
│   ├── variant_calling_bcftools.nf
│   ├── variant_calling_freebayes.nf
│   ├── variant_calling_gatk.nf
│   ├── annotation_vep.nf
│   ├── extract_cds_bed.nf
│   ├── extract_sequences.nf
│   ├── translate_to_protein.nf
│   └── identify_peptides.nf
├── conf/
│   └── nextflow.config            # Pipeline configuration
├── docker/
│   └── Dockerfile                 # Container definition for Docker profile
├── environment.yml                # Conda environment definition
└── README.md                      # Project documentation
</pre>
</details>

<details>
  <summary><b>Instalation</b></summary>

**Requirements**

- Nextflow ≥ 23.10
- Docker or Conda (depending of profile)
- Git

**Setup**

To set up the workflow, clone the repository and choose your preferred environment:
```bash
# Clone repository
git clone https://github.com/<your-username>/immunopeptidogenomics-workflow.git
cd immunopeptidogenomics-workflow

# Option 1: Use Docker
nextflow run main_SNV.nf -profile docker

# Option 2: Use Conda
nextflow run main_SNV.nf -profile conda
```
</pre>
</details>

<details>
  <summary><b>Input files</b></summary>

The workflow requires several input parameters that define sequencing data and reference files.

| **Parameter**        | **Description**                             | **Example**                                    |
|----------------------|---------------------------------------------|------------------------------------------------|
| `--reads`            | Path to FASTQ files (gzip)                  | `/<your-path>/SAMPLE_READS{1,2}.fastq.gz`   |
| `--ref_fasta`        | Reference genome in FASTA format            | `/<your-path>/GRCh38.fa`                    |
| `--prebuilt_bundle`  | Defines weather user can provide all refrence files along with corresponding index files, or index files need to be generated by pipeline|`true` or `false`|
| `--ref_fai`          | FASTA index                                 | `/<your-pathy>/GRCh38.fa.fai`                |
| `--ref_dict`         | Reference dictionary                        | `/<your-path>/GRCh38.dict`                  |
| `--hisat2_index`     | Path to HISAT2 index directory              | `/<your-path>/hisat2_index/`                |
| `--gff3_filtered`    | Annotation file in GFF3 formate (gzip)      | `/<your-path>/GRCh38.110.mainChr.gff3.gz`   |
| `--gff3_filtered_tbi`| Annotation file index                       | `/<your-path>/GRCh38.110.mainChr.gff3.gz.tbi`|
| `--cds_bed`          | Genomic coordinates for coding sequences    | `/<your-path>/cds.bed`                      |
| `--min_dp`           | Filtering value for the depth of reads coverage | `3`                                        |
| `--min_alt`          | Filtering value for the number of alternated nucleotides among the total number nucleotides covering the specific position | `3` |
| `--caller`           | Variant caller selection (`bcftools`, `freebayes`, `gatk`) | `bcftools`                      |
| `--outdir`           | Output directory                            | `./results_bcftools`                           |

If some of the parameters are not provided directly in the terminal, the pipeline will stop with an error message. 
</pre>
</details>

<details>
  <summary><b>Reference files</b></summary>

If user doesn’t already have the necessary reference genome and annotation files, they can be downloaded from publicly available sources.  
- Reference genome: [Ensembl GRCh38 Release 114](https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/) (Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)
- Gene model for  allignment and annotation: [Ensembl GTF](https://ftp.ensembl.org/pub/release-114/gff3/homo_sapiens/) (Homo_sapiens.GRCh38.114.chr.gff3.gz)

If .fai, .dict, or HISAT2 index files are missing, they will be automatically created by the pipeline module prep_reference.nf during execution.

Command-line Download Example:
```bash
# Create a directory for reference data and go their 
mkdir -p /<your-pathy>/genome && cd /<your-path>/genome

# Download reference genome: 
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Decompress file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download annotation file:
wget https://ftp.ensembl.org/pub/release-114/gff3/homo_sapiens/Homo_sapiens.GRCh38.114.chr.gff3.gz

# Decompress file
gunzip Homo_sapiens.GRCh38.114.gff3.gz

# Optional steps; requires pre-installation of specific tools (samtools, picard, hisat2, tabix)

# Generate fasta index (.fai):
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Generate sequence dictionary (.dict):
picard CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.primary_assembly.fa O=Homo_sapiens.GRCh38.dna.primary_assembly.dict

# Build HISAT2 genome index:
# Create output directory
mkdir -p hisat2_index
# Build index from reference genome
hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa hisat2_index/GRCh38

# Index the GFF3 Annotation
bgzip Homo_sapiens.GRCh38.114.chr.gff3
tabix -p gff Homo_sapiens.GRCh38.114.chr.gff3.gz
```
If you need to use another release, check the [Ensembl Archive List](https://asia.ensembl.org/Help/ArchiveList), which provides direct access to all historical Ensembl versions.

</pre>
</details>

<details>
  <summary><b>Output structure</b></summary>
<pre>
results/
├── alignment/
|   └── SampleName.sam
├── bam/ 
|   ├── SampleName.split.bam
|   ├── SampleName.fixmate.cs.bam  
|   ├── SampleName.markdup.bam
|   ├── SampleName.markdup.bam.bai 
|   ├── SampleName.sorted.bam   
│   └── SampleName.sorted.bam.bai
├── qc/fastp/
|   ├── SampleName.trim.fastq.gz
|   ├── SampleName.fastp.json  
│   └── SampleName.fastp.html            
└── vcf/
    ├── annotated_variants/  
    |      ├── SampleName.annotated.vcf.gz
    |      └── SampleName.annotated.vcf.gz.tbi
    ├── bcftools/
    |      ├── SampleName.raw.vcf.gz
    |      ├── SampleName.raw.vcf.gz.tbi
    |      ├── SampleName.filtered.vcf.gz
    |      └── SampleName.filtered.vcf.gz.tbi  
    ├── extract_cds_bed/ 
    |      ├── SampleName.affected_trascripts.txt
    |      ├── SampleName.cds_affected.bed
    |      ├── SampleName.protein_altering.bed
    |      ├── SampleName.protein_altering.positions.tsv
    |      ├── SampleName.protein_altering.vcf.gz
    |      └── SampleName.protein_altering.vcf.gz.tbi  
    ├── extract_sequences/ 
    |      ├── genome_mutated.fa
    |      ├── SampleName.mutated_cds.fa
    |      └── SampleName.original_cds.fa
    ├── protein_altering/ 
    |      ├── SampleName.affected_trascripts.txt
    |      ├── SampleName.protein_altering.bed
    |      ├── SampleName.protein_altering.positions.tsv
    |      ├── SampleName.protein_altering.vcf.gz
    |      └── SampleName.protein_altering.vcf.gz.tbi 
    └── translate_proteins/   
           ├── SampleName.mereged_cds.fa
           ├── SampleName.mutated_proteins.annot.txt
           └── SampleName.mutated_proteins.fa
                 
</pre>
</details>





