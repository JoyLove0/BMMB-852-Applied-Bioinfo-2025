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

### ------------------------ Variables ----------------------------
# How many reads to download N=ALL for all.
N ?= 1000000

# The number of CPUs to use.
NCPU ?= 4

# The URL for the Human transcriptome.
URL ?= http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# The name of the CDNA transcriptome file locally.
REF ?= refs/hsapiens/Homo_sapiens.GRCh38.cdna.all.fa.gz

# The name of the ensembl database to use for transcript to gene mapping.
ENSEMBL_DB ?= hsapiens_gene_ensembl

# The transcript to gene mapping file.
TX2GENE ?= refs/tx2gene_${ENSEMBL_DB}.csv

# The resulting counts file.
COUNTS ?= res/airway-counts.csv

# Final edger result
RESULT = res/edger.csv

# The design matrix.
DESIGN ?= airway.csv

# Check that required executables exist.
EXE = salmon fastq-dump parallel micromamba

# Flags we pass to parallel in each invocation
FLAGS = -v --eta --colsep , --header : -j ${NCPU}

# FASTQ read pair.
R1 ?= reads/read_1.fastq

# Second in pair.
R2 ?=

# Set MODE to SE if R2 is empty otherwise PE
ifeq ($(R2),)
	MODE ?= SE
else
	MODE ?= PE
endif

# The name of the sample that is processed
SAMPLE = $(basename $(notdir ${R1}))

# The name of the program to check
PROG = salmon

# Name of the output directory.
PREFIX ?= salmon/${SAMPLE}

# The name of the quantification file.
QUANT ?= ${PREFIX}/quant.sf

# The directory that holds the index.
IDX_DIR = $(dir ${REF})/idx

# The name of the index
IDX ?= ${IDX_DIR}/$(notdir ${REF}).salmon

# The index file name
IDX_FILE ?= ${IDX}/info.json

# Additional salmon options.
SALMON_FLAGS = -l A

### ------------------- Actual Workflow ----------------------------

# ------------------------------------------------------------
# Step 1: Get the data
# usage: make get_genome
# usage: make get_reads
# ------------------------------------------------------------
# Download and unzip the reference genome.
get_genome:
	mkdir -p $(dir ${REF})
	curl ${URL} > ${REF}
	seqkit stats refs/hsapiens/*.fa.gz
	
get_reads:
	bio search PRJNA229998 -H --csv --all > PRJNA229998.csv
	make -f src/workflows/airway.mk fastq
	seqkit stats reads/*.fastq

# ------------------------------------------------------------
# Step 2: Create the design file
# usage: make design
# ------------------------------------------------------------
design:
	make -f src/workflows/airway.mk design

# ------------------------------------------------------------
# Step 3: Index the genome
# usage: make index
# ------------------------------------------------------------
index:
	micromamba run -n salmon_env make -f src/run/salmon.mk index REF=${REF} NCPU=${NCPU}

# ------------------------------------------------------------
# Step 4: Align the reads
# usage: make align
# ------------------------------------------------------------
align: ${DESIGN}
	# Run the alignment on each sample.
	cat ${DESIGN} | parallel ${FLAGS} \
	                micromamba run -n salmon_env make -f src/run/salmon.mk run \
	                REF=${REF} \
	                SAMPLE={sample} \
	                R1=reads/{run}_1.fastq \
	                R2=reads/{run}_2.fastq

# ------------------------------------------------------------
# Step 5: Generate counts
# usage: make counts
# ------------------------------------------------------------
counts:
	make -f src/workflows/airway.mk counts

# ------------------------------------------------------------
# Step 6: Run edger analysis
# usage: make edger
# ------------------------------------------------------------
edger:
	make -f src/workflows/airway.mk edger
	cut -f 1,2,10 -d , res/edger.csv | head -11 > res/lowestFDR_expressed_report.txt
	
# ------------------------------------------------------------
# Step 7: Perform functional enrichment analysis
# usage: make enrich
# ------------------------------------------------------------
enrich:
	mkdir -p res/enrich
	# bio gprofiler -c res/edger.csv -d hsapiens
	cat gprofiler.csv | csvcut -c 2,3 | grep -iE 'vessel|lung|matrix' > res/enrich/biol_relav_report.txt
	mv gprofiler.csv res/enrich/gprofiler_full_results.csv
	@echo "# Generated functional enrichment report: biol_relav_report.txt"

### ------------------- End of Workflow ----------------------------
.PHONY: get_genome get_reads design index align counts edger enrich