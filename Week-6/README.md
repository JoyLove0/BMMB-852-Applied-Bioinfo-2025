# Week 6: Generate a BAM alignment file

## Make Makefile that includes rules for:
- Obtaining the genome
- Downloading sequencing reads from SRA
- index: Index the genome
- align: Generate a sorted and indexed BAM file by aligning reads to the genome
```
# Simple Makefile for genome alignment workflow

# Variables
GENOME_ACC = GCF_000882815.3
SRA_ACC = SRR9945583
REF = refs/ncbi_dataset/data/$(GENOME_ACC)/zika_paper_genome.fa
R1 = reads/$(SRA_ACC)_1.fastq
SAM = bam/$(SRA_ACC).sam
BAM = bam/$(SRA_ACC).bam

# Default target
all: align

# Step 1: Obtain the genome
genome:
	mkdir -p refs
	datasets download genome accession $(GENOME_ACC) --include genome,gff3,gtf --filename refs/zika_dataset.zip
	unzip -o refs/zika_dataset.zip -d refs
	mv refs/ncbi_dataset/data/$(GENOME_ACC)/*.fna $(REF)

# Step 2: Download reads from SRA
reads:
	mkdir -p reads
	fastq-dump -X 1500 --outdir reads --split-files $(SRA_ACC)

# Step 3: Index the genome
index: genome
	bwa index $(REF)

# Step 4: Align and generate sorted, indexed BAM file
align: index reads
	mkdir -p bam
	bwa mem $(REF) $(R1) > $(SAM)
	samtools sort $(SAM) > $(BAM)
	samtools index $(BAM)

# Clean up all generated files
clean:
	rm -rf refs reads bam
```
## Explain the use of the Makefile in your project.

### Overview of Zika Virus Genome Alignment Pipeline

This project automates a basic bioinformatics workflow using a **Makefile**. It downloads the Zika virus genome and sequencing reads, indexes the genome, and aligns the reads to produce a sorted and indexed BAM file.

---

### 📋 Workflow Overview

The pipeline includes the following steps:

1. **Obtain the reference genome** (Zika virus)
2. **Download sequencing reads** from the SRA
3. **Index the genome** using `bwa`
4. **Align reads** to the genome, sort, and index the output BAM

Each of these steps is defined as a separate **Makefile target**, allowing you to run them individually or all at once.

---

| **Target** | **Description**                                          |
| ---------- | -------------------------------------------------------- |
| `genome`   | Downloads and extracts the Zika virus genome             |
| `reads`    | Downloads reads from SRA                                 |
| `index`    | Indexes the reference genome using `bwa index`           |
| `align`    | Aligns reads, sorts SAM to BAM, and indexes the BAM file |
| `clean`    | Removes all generated files and folders                  |
| `all`      | Default target, runs everything up to `align`            |

### Usage
You can run the entire pipeline or individual steps using the make command.

### Run the Entire Pipeline
```
make
```
This is equivalent to:
```
make align
```
It will:
- Download the genome
- Download reads
- Index the genome
- Align reads and produce sorted + indexed BAM

### Step-by-Step Usage
1. ``` make genomen ```
Downloads and extracts the Zika virus genome.
Source: NCBI Genome via datasets CLI
Output: refs/zika_paper_genome.fa

2.``` make reads ```
Downloads sequencing reads from the SRA.
Accession: SRR9945583
Output: reads/SRR9945583_1.fastq
Only single-end reads are used in this example.

3. ```make index```
Indexes the reference genome using BWA.
Input: refs/zika_paper_genome.fa
Output: BWA index files (.bwt, .pac, etc.)
Note: make index depends on make genome.

4. ```make align```
Aligns reads to the genome, sorts the result, and creates an indexed BAM.
Input: Genome, reads
Output: bam/SRR9945583.bam and bam/SRR9945583.bam.bai
Note: make align depends on make index and make reads.

Extra: Clean Up Generated Files
``` make clean ```
Deletes all generated directories and files:
- refs/
- reads/
- bam/
  
### Customization
You can change the genome or read accession numbers by editing the top of the Makefile:
```
GENOME_ACC = GCF_000882815.3
SRA_ACC = SRR9945583
```
To process another sample, update SRA_ACC and re-run the relevant targets.

### File Structure After Running
.
├── Makefile
├── refs/
│   └── zika_paper_genome.fa
├── reads/
│   └── SRR9945583_1.fastq
└── bam/
    ├── SRR9945583.sam
    ├── SRR9945583.bam
    └── SRR9945583.bam.bai


## Visualize the resulting BAM files for both simulated reads and reads downloaded from SRA.
My SRA file aligned to nothing at first (see below). 
<img width="1440" height="900" alt="Screenshot 2025-10-05 at 11 17 43 PM" src="https://github.com/user-attachments/assets/38d5e4c7-ba5a-40d5-99ba-610209eecf4d" />

So, I had to download the whole SRA file for to get some hits. 

## Generate alignment statistics for the BAM file.
- What percentage of reads aligned to the genome?
- What was the expected average coverage?
- What is the observed average coverage?
- How much does the coverage vary across the genome? (Provide a visual estimate.)



