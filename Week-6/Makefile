# Makefile for my workflow

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
	echo "# Downloaded the Reference Genome with the accession number: $(GENOME_ACC). The path to the refernce is $(REF)."

# Step 2: Download reads from SRA
reads:
	mkdir -p reads
	fastq-dump -X 1500 --outdir reads --split-files $(SRA_ACC)
	echo "# Dowloaded Reads with the accession number: $(SRA_ACC)"

# Step 3: Index the genome
index: genome
	bwa index $(REF)
	echo "# Genome $(GENOME_ACC) was indexed."

# Step 4: Align and generate sorted, indexed BAM file
align: index reads
	mkdir -p bam
	bwa mem $(REF) $(R1) > $(SAM)
	samtools sort $(SAM) > $(BAM)
	samtools index $(BAM)
	echo "# Reads $(SRA_ACC) aligned to genome $(GENOME_ACC). The resulting sorted and indexed file is $(BAM)"

# Clean up all generated files
clean:
	rm -r refs reads bam
	echo " # ♫⋆｡♪ ₊˚♬ ﾟ.Clean Up! Clean Up! Everybody, Everywhere. Clean Up! Clean Up! Everybody do your share. ♫⋆｡♪ ₊˚♬ ﾟ.Directories refs, reads, and bam were deleted"
