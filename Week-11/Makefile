### Makefile for my workflow: analyzing Zika Virus data

### ------------------------ House Keeping ----------------------------
# Set the shell the commands run in.
SHELL = bash
# Execute all commands in a single shell.
.ONESHELL:
	# Run the shell with strict error checking.
	.SHELLFLAGS = -eu -o pipefail -c
# Delete target files if the command fails.
.DELETE_ON_ERROR:
# Warn if a variable is not defined.
MAKEFLAGS += --warn-undefined-variables
# Disable built-in rules.
MAKEFLAGS += --no-builtin-rules

### -------------------------- Variables ----------------------------
# All variables can be overridden from the command line.
GENOME_ACC ?= GCF_000882815.3
SRA_ACC ?= SRR9945583
REF ?= refs/zika_paper_genome.fa
R1 ?= reads/$(SRA_ACC)_1.fastq
R2 ?= reads/$(SRA_ACC)_2.fastq
SAM ?= bam/$(SRA_ACC).sam
BAM ?= bam/$(SRA_ACC).bam
VCF ?= vcf/variants.vcf.gz
GFF ?= refs/ncbi_dataset/data/$(GENOME_ACC)/zika_paper_genome.gff3
IGV_PORT ?= 60151

### -------------------- Actual Workflow ----------------------------
# Default target
all: genome reads_single reads_paired index align_single align_paired wiggle variants_sample variants_multi
	echo "# All steps completed for SRA: $(SRA_ACC) and Genome: $(GENOME_ACC)"

# ------------------------------------------------------------
# Step 1: Obtain the genome
# Usage: make genome GENOME_ACC=GCF_000882815.3
# ------------------------------------------------------------
genome:
	mkdir -p refs
	datasets download genome accession $(GENOME_ACC) --include genome,gff3,gtf --filename refs/zika_dataset.zip
	unzip -o refs/zika_dataset.zip -d refs
	mv refs/ncbi_dataset/data/$(GENOME_ACC)/*.fna refs/zika_paper_genome.fa
	echo "# Downloaded the Reference Genome with the accession number: $(GENOME_ACC). The path to the refernce is ${REF}"

# ------------------------------------------------------------
# Step 2: Download reads from SRA
# Step 2a: Single-end reads
# Usage: make reads_single SRR_ACC=SRR1972739 
# ------------------------------------------------------------
reads_single:
	mkdir -p reads
	fastq-dump -X 150000 --outdir reads --split-files $(SRA_ACC)
	fastqc reads/$(SRA_ACC)_1.fastq -o reads/
	echo "# Dowloaded Reads with the accession number: $(SRA_ACC)"
# ------------------------------------------------------------
# Step 2b: Paired-end reads
# Usage: make reads_paired SRR_ACC=SRR1972740
# ------------------------------------------------------------
reads_paired:
	mkdir -p reads
	fastq-dump -X 150000--outdir reads --split-files $(SRA_ACC)
	fastqc reads/$(SRA_ACC)_1.fastq reads/$(SRA_ACC)_2.fastq -o reads/
	echo "# Dowloaded Reads with the accession number: $(SRA_ACC)"

# ------------------------------------------------------------
# Step 3: Index the genome
# Usage: make index REF=ref/zika_paper_genome.fa
# ------------------------------------------------------------
index:
	bwa index $(REF)
	samtools faidx $(REF)
	echo "# Genome for $(GENOME_ACC) was indexed."

# ------------------------------------------------------------
# Step 4: Align and generate sorted, indexed BAM file
# Step 4a: Single-end reads
# Usage: make align_single SRR_ACC=SRR1972739 REF=ref/zika_paper_genome.fa
# ------------------------------------------------------------
align_single:
	mkdir -p bam
	bwa mem $(REF) reads/$(SRA_ACC)_1.fastq > bam/$(SRA_ACC).sam
	samtools sort bam/$(SRA_ACC).sam > bam/$(SRA_ACC).bam
	samtools index bam/$(SRA_ACC).bam
	samtools flagstat bam/$(SRA_ACC).bam > bam/$(SRA_ACC)_alignment_report.txt
	echo "# The resulting sorted and indexed file is bam/$(SRA_ACC).bam"
# ------------------------------------------------------------
# Step 4b: Paired-end reads
# Usage: make align_paired SRR_ACC=SRR1972740 REF=ref/zika_paper_genome.fa
# ------------------------------------------------------------
align_paired:
	mkdir -p bam
	bwa mem $(REF) reads/$(SRA_ACC)_1.fastq reads/$(SRA_ACC)_2.fastq > bam/$(SRA_ACC).sam
	samtools sort bam/$(SRA_ACC).sam > bam/$(SRA_ACC).bam
	samtools index bam/$(SRA_ACC).bam
	samtools flagstat bam/$(SRA_ACC).bam > bam/$(SRA_ACC)_alignment_report.txt
	echo "# The resulting sorted and indexed file is bam/$(SRA_ACC).bam"

# ------------------------------------------------------------
# Step 6: Variant calling + automated visualization for one sample
# Usage: make variants SRR_ACC=SRR1972739 REF=ref/zika_paper_genome.fa  F1='-d 100 --annotate ...' F2='--ploidy 2 --annotate ...' IGV_PORT=60151
# ------------------------------------------------------------
variants_sample:
	mkdir -p vcf
	echo "# Calling variants for sample $(SRA_ACC)..."
	bcftools mpileup -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP' -O u -f ${REF} bam/$(SRA_ACC).bam | \
	bcftools call --ploidy 2 --annotate 'FORMAT/GQ' -mv -O u | \
	bcftools norm -f ${REF} -d all -O u | \
	bcftools sort -O z > vcf/$(SRA_ACC)_variants.vcf.gz

# ------------------------------------------------------------
# optional clean up command
# Usage: make clean
# ------------------------------------------------------------
clean:
	rm -r refs reads bam vcf
	echo " # ♫⋆｡♪ ₊˚♬ ﾟ.Clean Up! Clean Up! Everybody, Everywhere. Clean Up! Clean Up! Everybody do your share. ♫⋆｡♪ ₊˚♬ ﾟ.Directories refs, reads, and bam were deleted"