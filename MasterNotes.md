# Metagenomic Analysis Pipeline


## Dataset
- Source: European Nucleotide Archive (ENA), European Bioinformatics Institue (EMBL-EBI), 2017.
- Study Accession: PRJEB21528
- Control Sample 01:
  - Sample Accession - SAMEA104142370
  - Run Accession - ERR2017415
- Control Sample 02:
  - Sample Accession - SAMEA104142365
  - Run Accession - ERR2017420.
- ACVD Sample 01:
  - Sample Accession - SAMEA104142287
  - Run Accession - ERR2017598.
- ACVD Sample 02:
  - Sample Accession - SAMEA104142288
  - Run Accession - ERR2017599.

## Workflow Overview
The following workflow was followed for all 4 samples, and the sequence of commands was consistent. For the sake of conserving space, the code for the first control is reproduced below.
1. Download raw sequencing data
2. Assess raw read quality with FastQC
3. Trim adapters and low quality bases with Trimmomatic
4. Assess trimmed read quality with FastQC
5. Identify taxonomies present in our dataset using Kraken2
6. Visualize diversity and abundance using R

## Step 1: Set up project directories and download sequencing reads

The first step is to create a reproducible project structure so that raw reads, trimmed reads, logs, and outputs are separated. This makes the workflow easier to follow and reduces path errors in downstream scripts. We then download the sequencing data from the ENA database.

-------------------
```bash
# Create a main project directory to store all files for this workflow
mkdir final_project

# Move into the project directory so all subsequent files are organized here
cd final_project

# Create a subdirectory to store adapter sequences used during trimming
mkdir adapters

# Create a subdirectory for sequencing reads and FastQC reports
mkdir raw_reads

# Create a subdirectory for trimmed reads and FastQC reports
mkdir trimmed_reads

# Create a subdirectory to store log files
mkdir logs

# Create a subdirectory to write slurm scripts
mkdir scripts

# Create a subdirectory to store a database for Kraken2
mkdir kraken_db

# Create a subdirectory to store the outputs from Kraken2
mkdir kraken_output

# Download raw read files from publicly available dataset
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR201/005/ERR2017415/ERR2017415_1.fastq.gz
$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR201/005/ERR2017415/ERR2017415_2.fastq.gz
