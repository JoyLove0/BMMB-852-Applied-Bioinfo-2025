### Makefile for my workflow: analyzing Zika Virus data
# The result should consist of:
# - A genome named in a user-friendly name
# - FASTQ read data named by the samples
# - FASTQC reports for each read
# - Alignments and coverage files in BAM and BW formats.
# - A statistics alignment report for BAM file
# - A VCF file with variants called from the BAM file

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
REF ?= refs/ncbi_dataset/data/$(GENOME_ACC)/zika_paper_genome.fa
R1 ?= reads/$(SRA_ACC)_1.fastq
R2 ?= reads/$(SRA_ACC)_2.fastq
SAM ?= bam/$(SRA_ACC).sam
BAM ?= bam/$(SRA_ACC).bam
VCF ?= vcf/variants.vcf.gz
GFF ?= refs/ncbi_dataset/data/$(GENOME_ACC)/zika_paper_genome.gff3
F1 ?= -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP' # Flags for the mpileup command
F2 ?= --ploidy 2 --annotate 'FORMAT/GQ' # Flags for the call command
IGV_PORT ?= 60151 # IGV network interface port (default: 60151)

### -------------------- Actual Workflow ----------------------------
# Default target
all: genome reads index align wiggle
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
	echo "# Downloaded the Reference Genome with the accession number: $(GENOME_ACC). The path to the refernce is refs/zika_paper_genome.fa."

# ------------------------------------------------------------
# Step 2: Download reads from SRA
# Step 2a: Single-end reads
# Usage: make reads_single SRR_ACC=SRR1972739 
# ------------------------------------------------------------
reads_single:
	mkdir -p reads
	fastq-dump -X 1500 --outdir reads --split-files $(SRA_ACC)
	fastqc $(SRA_ACC)_1.fastq -o reads/
	echo "# Dowloaded Reads with the accession number: $(SRA_ACC)"
# ------------------------------------------------------------
# Step 2b: Paired-end reads
# Usage: make reads_paired SRR_ACC=SRR1972740
# ------------------------------------------------------------
reads_paired:
	mkdir -p reads
	fastq-dump -X 1500 --outdir reads --split-files $(SRA_ACC)
	fastqc reads/$(SRA_ACC)_1.fastq reads/$(SRA_ACC)_2.fastq -o reads/
	echo "# Dowloaded Reads with the accession number: $(SRA_ACC)"

# ------------------------------------------------------------
# Step 3: Index the genome
# Usage: make index REF=ref/zika_paper_genome.fa
# ------------------------------------------------------------
index:
	bwa index $(REF)
	echo "# Genome $(GENOME_ACC) was indexed."

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
# Step 5: Coverts the BAM file into a BigWig file
# Usage: make wiggle SRR_ACC=SRR1972739 REF=ref/zika_paper_genome.fa
# ------------------------------------------------------------
wiggle:
	samtools faidx ${REF}
	LC_ALL=C; bedtools genomecov -ibam  bam/$(SRA_ACC).bam -split -bg | \
		sort -k1,1 -k2,2n > bam/$(SRA_ACC).bedgraph
	bedGraphToBigWig bam/$(SRA_ACC).bedgraph ${REF} bam/$(SRA_ACC).bw
	echo "# The BAM file bam/$(SRA_ACC).bam was converted to a BigWig file bam/$(SRA_ACC).bw"

# ------------------------------------------------------------
# Step 6: Variant calling + automated visualization for one sample
# Usage: make variants SRR_ACC=SRR1972739 REF=ref/zika_paper_genome.fa  F1='-d 100 --annotate ...' F2='--ploidy 2 --annotate ...' IGV_PORT=60151
# ------------------------------------------------------------
variants_sample:
	mkdir -p vcf
	echo "# Calling variants for sample $(SRA_ACC)..."
	bcftools mpileup ${F1} -O u -f ${REF} bam/$(SRA_ACC).bam | \
		bcftools call ${F2} -m -v -O z -o vcf/$(SRA_ACC).vcf.gz
	bcftools index vcf/$(SRA_ACC).vcf.gz
	echo "# VCF file created and indexed: vcf/$(SRA_ACC).vcf.gz"
	echo ""
	echo "# Visualization via IGV network interface (port $(IGV_PORT))"
	echo "# Make sure IGV is running with port $(IGV_PORT) enabled."
	echo "new" | nc localhost $(IGV_PORT) || true
	echo "genome $(REF)" | nc localhost $(IGV_PORT) || true
	echo "load $(GFF)" | nc localhost $(IGV_PORT) || true
	echo "load bam/$(SRA_ACC).bam" | nc localhost $(IGV_PORT) || true
	echo "load vcf/$(SRA_ACC).vcf.gz" | nc localhost $(IGV_PORT) || true
	echo "# IGV session loaded for sample $(SRA_ACC)."

# ------------------------------------------------------------
# Step 7: Multisample VCF creation + automated visualization
# Usage: make variants_multi IGV_PORT=60151 VCF_LIST='vcf/SRR1.vcf.gz vcf/SRR2.vcf.gz'
# ------------------------------------------------------------
variants_multi:
	mkdir -p vcf
	echo "# Merging VCF files into a multisample VCF..."
	bcftools merge -Oz -o vcf/multisample.vcf.gz $(VCF_LIST)
	bcftools index vcf/multisample.vcf.gz
	echo "# Multisample VCF created: vcf/multisample.vcf.gz"
	echo ""
	echo "# Visualization via IGV network interface (port $(IGV_PORT))"
	echo "# Make sure IGV is running with port $(IGV_PORT) enabled."
	echo "new" | nc localhost $(IGV_PORT) || true
	echo "genome $(REF)" | nc localhost $(IGV_PORT) || true
	echo "load $(GFF)" | nc localhost $(IGV_PORT) || true
	for BAM in bam/*.bam; do \
		echo "load $$BAM" | nc localhost $(IGV_PORT) || true; \
	done
	echo "load vcf/multisample.vcf.gz" | nc localhost $(IGV_PORT) || true
	echo "# ✅ IGV session loaded with multisample data."

# ------------------------------------------------------------
# optional clean up command
# Usage: make clean
# ------------------------------------------------------------
clean:
	rm -r refs reads bam
	echo " # ♫⋆｡♪ ₊˚♬ ﾟ.Clean Up! Clean Up! Everybody, Everywhere. Clean Up! Clean Up! Everybody do your share. ♫⋆｡♪ ₊˚♬ ﾟ.Directories refs, reads, and bam were deleted"