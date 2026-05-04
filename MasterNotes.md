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

### Slurm Script: 
Parameter choices (ILLUMINACLIP, SLIDINGWINDOW, MINLEN) were selected to balance read quality and data retention. Adapter clipping removes non-biological sequences, SLIDINGWINDOW:4:20 trims low-quality regions based on a common Q20 threshold (≈1% error rate) in a window of 4 bases, and MINLEN:50 ensures reads are long enough (>50 bases) to be useful for assembly while discarding overly short, unreliable fragments. The script is produced below.

-------------------
```bash
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

The slurm script was submitted as a job to the HPC scheduler to run on a compute node.

--------------------
```bash
sbatch trim_reads.slurm
```
--------------------

80.16% of reads survived in both forward and reverse trimming.

--------------------

## Step 4: Assess trimmed read quality with FastQC

-------------------
```bash
# Move back into the main project directory and then into the subdirectory with the trimmed reads
cd ..
cd trimmed_reads

# Run FastQC on paired trimmed reads
srun --pty bash
module load fastqc
fastqc *_paired.fastq.gz
```
-------------------

## Step 5: Identify taxonomies present in our dataset using Kraken2
Kraken2 uses trimmed reads to perform k-mer-based taxonomic classification. Kraken2 requires a databse, and we downloaded a publicly available standard database.

-------------------
```bash
# Create an environment to activate kraken2
module load anaconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda create -n kraken2 -c bioconda kraken2 -y
conda activate kraken2

# Download a standard database for kraken2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260226.tar.gz
tar -xzvf k2_standard_08_GB_20260226.tar.gz

# Run kraken2 on paired trimmed reads
kraken2 \
--db /home/mjd346/final_project/kraken/ \
--threads 8 \
--paired \
--output /home/mjd356/final_project/kraken_output/paired.output \
--report /home/mjd356/final_project/kraken_output/paired.report \
--use-names \
/home/mjd356/final_project/trimmed_reads/ERR2017420_forward_paired.fastq.gz \
/home/mjd356/final_project/trimmed_reads/ERR2017420_reverse_paired.fastq.gz
```
-------------------

## Step 6: Visualize diversity and abundance using R
We created a relative abundance graph and a UMAP in R to visualize the diversity and abundance of microbes in our samples.

Public links to R workflows:
- Relative Abundance Graph from Presentation: https://rpubs.com/skar/1425504
- UMAP: https://rpubs.com/skar/1425493
- Relative Abundance Graph code used in writeup addressing changes.

#install.packages(c("dplyr", "tidyr", "ggplot2"))

#library(dplyr)
#library(tidyr)
#library(ggplot2)

#file names
files <- c(
  "control_01.report",
  "control_02.report",
  "treatment_01.report",
  "treatment_02.report"   # fix name if needed
)

#read through kraken files
parse_kraken <- function(file) {
  lines <- readLines(file)
  
  parsed <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    
    if (length(parts) < 6) return(NULL)
    
    data.frame(
      percent = as.numeric(parts[1]),
      rank = parts[4],
      name = trimws(parts[6]),
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, parsed)
  
  # Keep genus only
  df <- df[df$rank == "G", c("name", "percent")]
  colnames(df) <- c("Taxon", "Percent")
  
  return(df)
}

#read files

data_list <- lapply(seq_along(files), function(i) {
  df <- parse_kraken(files[i])
  colnames(df)[2] <- files[i]
  df
})

# merge samples

merged <- Reduce(function(x, y) merge(x, y, by = "Taxon", all = TRUE), data_list)

# Replace NA with 0
merged[is.na(merged)] <- 0

#matrix
merged_matrix <- merged
rownames(merged_matrix) <- merged_matrix$Taxon
merged_matrix$Taxon <- NULL

#convert into long format
full_df <- merged_matrix %>%
  as.data.frame() %>%
  mutate(Taxon = rownames(.)) %>%
  pivot_longer(
    cols = -Taxon,
    names_to = "sample",
    values_to = "percent"
  )

#normalize data

full_df <- full_df %>%
  group_by(sample) %>%
  mutate(percent = percent / sum(percent) * 100) %>%
  ungroup()

#collapse 1% into one another for an other category
full_df <- full_df %>%
  mutate(
    Taxon = ifelse(percent < 1, "Other", Taxon)
  )


plot_df <- full_df %>%
  group_by(sample, Taxon) %>%
  summarise(percent = sum(percent), .groups = "drop")

#order by importance
plot_df <- plot_df %>%
  group_by(Taxon) %>%
  mutate(total = sum(percent)) %>%
  ungroup() %>%
  arrange(desc(total))

plot_df$Taxon <- factor(plot_df$Taxon, levels = unique(plot_df$Taxon))

#plot
ggplot(plot_df, aes(x = sample, y = percent, fill = Taxon)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Relative Abundance of Bacterial Genera (All Taxa, <1% Grouped)",
    x = "Sample",
    y = "Relative Abundance (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )  
