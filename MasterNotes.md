# Metagenomic Analysis Pipeline


## Dataset
- Source: European Nucleotide Archive (ENA), European Bioinformatics Institue (EMBL-EBI), 2017.
- Study Accession: PRJEB21528
- Control Sample 01:
  - Sample Accession - SAMEA104142370
  - Run Accession - ERR2017415
- Control Sample 02:
  - Sample Accession - SAMEA104142365
  - Run Accession - ERR2017420
- ACVD Sample 01:
  - Sample Accession - SAMEA104142287
  - Run Accession - ERR2017598
- ACVD Sample 02:
  - Sample Accession - SAMEA104142288
  - Run Accession - ERR2017599

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

# Move into the subdirectory to store raw reads
cd raw_reads

# Download raw read files from publicly available dataset
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR201/005/ERR2017415/ERR2017415_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR201/005/ERR2017415/ERR2017415_2.fastq.gz
```
--------------------

## Step 2: Assess raw read quality with FastQC

FastQC is used to evaluate the quality of the raw sequencing reads before any downstream analysis. This is important because poor-quality bases, adapter contamination, and abnormal sequence duplication can negatively affect trimming, assembly, and viral detection. By inspecting the raw reads first, we can justify why trimming is needed.

-------------------
```bash
srun --pty bash
module load fastqc
fastqc ERR2017415_*
```
-------------------

## Step 3: Trim adapters and low-quality bases with Trimmomatic

Trimming removes adapter contamination and low-quality sequence from the reads. We use paired-end trimming so that the relationship between the two reads from the same DNA fragment is preserved whenever possible.

-------------------
```bash
# Move back into the main project directory and then into the subdirectory for scripts
cd ..
cd scripts

# Open a text file for the script
nano trim_control_01
```
-------------------

### Slurm Script: Parameter choices (ILLUMINACLIP, SLIDINGWINDOW, MINLEN) were selected to balance read quality and data retention. Adapter clipping removes non-biological sequences, SLIDINGWINDOW:4:20 trims low-quality regions based on a common Q20 threshold (≈1% error rate) in a window of 4 bases, and MINLEN:50 ensures reads are long enough (>50 bases) to be useful for assembly while discarding overly short, unreliable fragments.

-------------------
```bash
# The script is produced below
#!/bin/bash
#SBATCH --job-name=trim_control_01
#SBATCH --output=/home/mjd356/final_project/logs/trim_%j.out
#SBATCH --error=/home/mjd356/final_project/logs/trim_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjd356@georgetown.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G

shopt -s expand_aliases
module load trimmomatic

BASE=/home/mjd356/final_project
R1=$BASE/raw_reads/ERR2017415_1.fastq.gz
R2=$BASE/raw_reads/ERR2017415_2.fastq.gz
ADAPTERS=$BASE/adapters/TruSeq3-PE.fa
OUTDIR=$BASE/trimmed_reads

trimmomatic PE -threads ${SLURM_CPUS_PER_TASK} \
  $R1 $R2 \
  $OUTDIR/ERR2017415_forward_paired.fastq.gz $OUTDIR/ERR2017415_forward_unpaired.fastq.gz \
  $OUTDIR/ERR2017415_reverse_paired.fastq.gz $OUTDIR/ERR2017415_reverse_unpaired.fastq.gz \
  ILLUMINACLIP:$ADAPTERS:2:30:10 \
  SLIDINGWINDOW:4:20 \
  MINLEN:50
```
-------------------
