### Makefile for my workflow

### ---- House Keeping ----
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

### -----Actual Workflow-----
# Default target
all: genome reads index align wiggle
	echo "# All steps completed for SRA: $(SRA_ACC) and Genome: $(GENOME_ACC)"

# Step 1: Obtain the genome
genome:
	mkdir -p refs
	datasets download genome accession $(GENOME_ACC) --include genome,gff3,gtf --filename refs/zika_dataset.zip
	unzip -o refs/zika_dataset.zip -d refs
	mv refs/ncbi_dataset/data/$(GENOME_ACC)/*.fna refs/zika_paper_genome.fa
	echo "# Downloaded the Reference Genome with the accession number: $(GENOME_ACC). The path to the refernce is $(REF)."

# Step 2: Download reads from SRA
reads:
	mkdir -p reads
	fastq-dump -X 1500 --outdir reads --split-files $(SRA_ACC)
	echo "# Dowloaded Reads with the accession number: $(SRA_ACC)"

# Step 3: Index the genome
index:
	bwa index $(REF)
	echo "# Genome $(GENOME_ACC) was indexed."

# Step 4: Align and generate sorted, indexed BAM file
align:
	mkdir -p bam
	bwa mem $(REF) reads/$(SRA_ACC)_1.fastq > bam/$(SRA_ACC).sam
	samtools sort bam/$(SRA_ACC).sam > bam/$(SRA_ACC).bam
	samtools index bam/$(SRA_ACC).bam
	samtools flagstat bam/$(SRA_ACC).bam > bam/$(SRA_ACC)_alignment_report.txt
	echo "# Reads $(SRA_ACC) aligned to genome $(GENOME_ACC). The resulting sorted and indexed file is bam/$(SRA_ACC).bam"

# Step 5: Coverts the BAM file into a BigWig file
wiggle:
	samtools faidx ${REF}
	LC_ALL=C; bedtools genomecov -ibam  bam/$(SRA_ACC).bam -split -bg | \
		sort -k1,1 -k2,2n > bam/$(SRA_ACC).bedgraph
	bedGraphToBigWig bam/$(SRA_ACC).bedgraph ${REF} bam/$(SRA_ACC).bw
	echo "# The BAM file bam/$(SRA_ACC).bam was converted to a BigWig file bam/$(SRA_ACC).bw"

clean:
	rm -r refs reads bam
	echo " # ♫⋆｡♪ ₊˚♬ ﾟ.Clean Up! Clean Up! Everybody, Everywhere. Clean Up! Clean Up! Everybody do your share. ♫⋆｡♪ ₊˚♬ ﾟ.Directories refs, reads, and bam were deleted"
