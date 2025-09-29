# Week 5: Obtain and visualize FASTQ data from SRA

## Write a Bash script:
- Reuse and expand your code from last week.
- Create a bash shell script with the code from last week.
- Add commands to download at least one sequencing dataset using the SRR number(s).
- Download only a subset of the data that would provide approximately 10x genome coverage. Briefly explain how you estimated the amount of data needed for 10x coverage.

```
: '
Goals: The objective of this code is to download sequence
data from a Short Read Archive SRR accession number. This 
code is intended to reproducible via user plug-ins. 
'
### ====================== MAKE CODE RUN SMOOTHER ====================== ###
set -x

### ====================== USER PLUG-INS ====================== ###
genbank_accession_number="GCF_000882815.3"  # Zika genome
srr_accession_number="SRR3194431"           # SRR ID
wd="reads"                             # Output folder for reads

### ================ DEPENDENCIES (Run once) ================= ###
# micromamba activate bioinfo

### Dependencies (do this once)
# micromamba activate bioinfo

###################################### NO EDIT BEYOND THIS LINE ################################
### Last Weeks Code:  Download genome (FASTA + annotation)
#datasets download genome accession $genebank_accession_number --include genome,gff3,gtf
#unzip ncbi_dataset.zip
#cd ncbi_dataset/data/GCF_000882815.3/
#mv GCF_000882815.3_ViralProj36615_genomic.fna zika_paper_genome.fa

### New Code: Create directory and download 1500 reads
mkdir -p $wd 
fastq-dump -X 1500 -F --outdir $wd --split-files $srr_accession_number
cd $wd

### ========== COVERAGE CALCULATION ========== ###
: '
Coverage = Total bases sequenced / Genome size

Goal: 10× coverage
Zika virus genome size = 10,794 bp
Total bases needed = 10 × 10,794 = 107,940 bases

From SRR3194431 metadata (Command: bio search SRR3194431):
- Total sequenced bases = 5,008,831,278
- Total reads = 66,528,035
- Average read length ≈ 75 bp

Therefore, full dataset provides:
Coverage = 5,008,831,278 / 10,794 ≈ 464,041×

To achieve only 10× coverage, we need:
107,940 bases ÷ 75 bp ≈ 1,439 reads

This script downloads a small subset (1,000–1,500 reads) for test purposes,
which is roughly enough for 7–10× coverage, assuming uniform distribution and no host contamination.

NOTE: Since this is RNA-Seq from human cells infected with Zika virus,
not all reads will be from the viral genome.
Actual effective coverage of the Zika genome may be lower.
'

### Note that to run this script use: bash ssr_download.sh
```

## Quality assessment:
- Generate basic statistics on the downloaded reads (e.g., number of reads, total bases, average read length).
- Run FASTQC on the downloaded data to generate a quality report.
- Evaluate the FASTQC report and summarize your findings.

## Compare sequencing platforms:
- Search the SRA for another dataset for the same genome, but generated using a different sequencing platform (e.g., if original data was Illumina select PacBio or Oxford Nanopore).
- Briefly compare the quality or characteristics of the datasets from the two platforms.
