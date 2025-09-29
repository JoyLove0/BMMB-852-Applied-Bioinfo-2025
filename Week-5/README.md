# Week 5: Obtain and visualize FASTQ data from SRA

## Write a Bash script:
- Reuse and expand your code from last week.
- Create a bash shell script with the code from last week.
- Add commands to download at least one sequencing dataset using the SRR number(s).
- Download only a subset of the data that would provide approximately 10x genome coverage. Briefly explain how you estimated the amount of data needed for 10x coverage.

Code:
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
Command:
```
bash ssr_download.sh 
```
Ouput:
```
+ genbank_accession_number=GCF_000882815.3
+ srr_accession_number=SRR3194431
+ wd=reads
+ mkdir -p reads
+ fastq-dump -X 1500 -F --outdir reads --split-files SRR3194431
Rejected 1500 READS because READLEN < 1
Read 1500 spots for SRR3194431
Written 1500 spots for SRR3194431
+ cd reads
+ : '
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
```
## Quality assessment:
- Generate basic statistics on the downloaded reads (e.g., number of reads, total bases, average read length).

Command:
```
seqkit stats SRR3194431_1.fastq
```

Output:
```
file                format  type  num_seqs  sum_len  min_len  avg_len  max_len
SRR3194431_1.fastq  FASTQ   DNA      1,500  113,166       35     75.4       76
```
- Run FASTQC on the downloaded data to generate a quality report.

Command:
```
fastqc SRR3194431_1.fastq
```

Output:
```
null
Started analysis of SRR3194431_1.fastq
Approx 65% complete for SRR3194431_1.fastq
Analysis complete for SRR3194431_1.fastq
```
- Evaluate the FASTQC report and summarize your findings.
<img width="1440" height="900" alt="Screenshot 2025-09-28 at 10 43 28 PM" src="https://github.com/user-attachments/assets/41d1c567-a16c-4b10-9af9-fcd1f77ad9ca" />

<img width="1440" height="900" alt="Screenshot 2025-09-28 at 10 43 28 PM" src="https://github.com/user-attachments/assets/cdef2faf-f1b8-4ba4-b996-cb90fa51516e" />

Overall, the data looks like we can do some really fun analysis after some quality control. There is a signficant portion of the data is in the "green", meaning they have a phred score between 30 to 36. This give me more condidence in our base calls as a q-score of 30 means there is a 1 in 1,000 chance of an incorrect base call (99.9% accuracy). I would still, however, reccomend to trim this data using something like Trimmomatic because there is significant portion of the data still in the red as indicated by the whiskers of the box plot.

## Compare sequencing platforms:
- Search the SRA for another dataset for the same genome, but generated using a different sequencing platform (e.g., if original data was Illumina select PacBio or Oxford Nanopore).
- Briefly compare the quality or characteristics of the datasets from the two platforms.

While navigating ncbi using this [tool](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=1&WebEnv=MCID_68d9daa17e6bb94164f643c1&o=instrument_s%3Ad), the following Zika SRA stood out: SRR9945583. It was sequenced by a next-generation sequencing (NGS) platform developed by BGI (now known as MGI) that competes with Illumina's platforms. 

More info [here](https://www.ncbi.nlm.nih.gov/sra/SRX6694077[accn]).

I used the same proccess for this data as above:
Command
```
bash ssr_download.sh # srr_accession_number="SRR9945583" now
seqkit stats SRR9945583_1.fastq
fastqc SRR9945583_1.fastq 
```

Outputs:
```
file                format  type  num_seqs  sum_len  min_len  avg_len  max_len
SRR9945583_1.fastq  FASTQ   DNA      1,500   37,381       18     24.9       42
```
<img width="551" height="324" alt="Screenshot 2025-09-28 at 10 46 00 PM" src="https://github.com/user-attachments/assets/a67da09a-2797-4880-9186-c95f6f294e37" />

<img width="2880" height="1800" alt="Screenshot 2025-09-28 at 11 33 25 PM" src="https://github.com/user-attachments/assets/eaa2200f-00e1-4c87-a926-ea6c591843ed" />

NOTE: General differences between sequencing technology:

| Feature                  | **BGISEQ-500**                                  | **Illumina NextSeq 500**                       |
| ------------------------ | ----------------------------------------------- | ---------------------------------------------- |
| **Manufacturer**         | BGI (Beijing Genomics Institute)                | Illumina (USA)                                 |
| **Sequencing Chemistry** | cPAS (combinatorial Probe-Anchor Synthesis)     | SBS (Sequencing by Synthesis)                  |
| **Read Type**            | Single-end or paired-end                        | Single-end or paired-end                       |
| **Typical Read Lengths** | 50 bp, 100 bp, 150 bp                           | 75 bp, 150 bp                                  |
| **Output per Run**       | Up to ~200 Gb                                   | Up to ~120 Gb                                  |
| **Flow Cell Design**     | Patterned array with DNA nanoballs (DNB)        | Flow cell clusters                             |
| **Error Profile**        | Slightly higher indel rate (vs Illumina)        | Lower error rate, especially for substitutions |
| **Instrument Size/Cost** | Lower instrument and reagent costs              | Higher cost overall                            |
| **Adoption / Ecosystem** | Less common outside China                       | Industry standard globally                     |
| **Data Compatibility**   | Mostly compatible with Illumina-style pipelines | Fully supported across most tools              |

The BGISEQ-500 dataset also starts off strong, with many reads in the high-quality zone. However, there's a noticeable quality drop around 30 bp, and the average read length is significantly shorter (~42 bp) due to trimming or lower-quality sequencing toward the 3' end. This suggests the need for more aggressive trimming and possibly read length filtering during preprocessing. While both datasets are usable, the Illumina data offers more consistent quality across the full read length, which may lead to more confident downstream analysis — especially in applications like variant calling or transcript quantification. See the table below for summary of the datasets:

| Feature                        | **Illumina NextSeq 500**                                                             | **BGISEQ-500**                                                              |
| ------------------------------ | ------------------------------------------------------------------------------------ | --------------------------------------------------------------------------- |
| **Overall Quality**            | Generally high — large portion of reads have **Phred scores 30–36**                  | Also strong early on, with a good portion in the **green zone**             |
| **Phred Score Interpretation** | Q30 = 99.9% accuracy → very reliable base calling                                    | Same scale used — Q30 still indicates high-confidence calls                 |
| **Notable Quality Features**   | - Long stretches of high-quality bases  <br> - Whiskers indicate lower-quality tails | - Sharp **quality drop around 30 bp** <br> - Lower quality after that point |
| **Average Read Length**        | ~75 bp (or user-defined during sequencing)                                           | Dropped to **~42 bp** due to quality-based truncation                       |
| **Recommended QC Steps**       | **Trimming** (e.g., with *Trimmomatic*) to remove low-quality reads                  |  **Aggressive trimming** advised, especially after 30 bp                   |
| **Confidence in Data**         | High — consistent with expected Illumina quality                                     | Moderate — early bases are good, but **quality decay is steeper**           |

