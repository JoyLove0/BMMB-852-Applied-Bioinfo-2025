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

REF          ?= refs/chr22.genome.fa
GTF          ?= refs/chr22.gtf
READS        ?= reads
BAM      ?= bam
DESIGN       ?= design.csv
FLAGS        ?= -v --eta --colsep , --header : -j 4
COUNTS       ?= counts_final.txt
COUNTS_TXT   ?= counts.txt
STR_COUNTS   ?= counts_strand_specific.txt
FORM_COUNTS  ?= counts_formated.txt
TX2GENE      ?= tx2gene.tsv

R1 ?= ${READS}/{sample}_1.fq

### ------------------- Actual Workflow ----------------------------

# ------------------------------------------------------------
# Step 1: Get the data
# usage: make get_data
# ------------------------------------------------------------
get_data:
	mkdir -p refs reads
	wget -nc http://data.biostarhandbook.com/data/uhr-hbr.tar.gz
	tar xzvf uhr-hbr.tar.gz
	seqkit stats refs/*.fa
	seqkit stats reads/*.fq

# ------------------------------------------------------------
# Step 2: Create the design file
# usage: make design
# ------------------------------------------------------------
design:
	@cat << EOF > ${DESIGN}
	sample,group
	HBR_1,HBR
	HBR_2,HBR
	HBR_3,HBR
	UHR_1,UHR
	UHR_2,UHR
	UHR_3,UHR
	EOF
	@echo "# Created design file: ${DESIGN}"

# ------------------------------------------------------------
# Step 3: Index the genome
# usage: make index
# ------------------------------------------------------------
index:
	make -f src/run/hisat2.mk index REF=${REF}

# ------------------------------------------------------------
# Step 4: Align the reads
# usage: make align
# ------------------------------------------------------------
align: ${DESIGN}
	mkdir -p ${BAM}
	cat ${DESIGN} | parallel ${FLAGS} \
		make -f src/run/hisat2.mk REF=${REF} SAMPLE={sample} R1=${READS}/{sample}_R1.fq BAM=${BAM}/{sample}.bam run

# ------------------------------------------------------------
# Step 5: Generate counts
# usage: make counts
# ------------------------------------------------------------
counts: $(COUNTS_TXT) $(STR_COUNTS) $(FORM_COUNTS)

${COUNTS_TXT}:
	featureCounts -a ${GTF} -o ${COUNTS_TXT} \
		${BAM}/*.bam

${STR_COUNTS}:
	featureCounts -a ${GTF} -o ${STR_COUNTS} -s 2 \
		${BAM}/*.bam

${FORM_COUNTS}: ${COUNTS_TXT}
	Rscript src/r/format_featurecounts.r -c ${COUNTS_TXT} -o ${FORM_COUNTS}

### ---------------- TX2GENE creation ----------------

${TX2GENE}:
	Rscript src/r/create_tx2gene.r -d hsapiens_gene_ensembl -o ${TX2GENE}

### ---------------- Combine final counts ----------------

${COUNTS}: ${DESIGN} ${TX2GENE} ${FORM_COUNTS}
	micromamba run -n stats Rscript src/r/combine_featurecounts.r \
		-c ${FORM_COUNTS} -t ${TX2GENE} -o ${COUNTS}

# ------------------------------------------------------------