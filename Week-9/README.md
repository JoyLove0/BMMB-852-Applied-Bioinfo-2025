# Week 9: Revising and Improving Your Automation Code
---
# Here are the revisions:
- There is a typo (a missing last number) in bio search PRJNA31329 that's why it does not return values
- Just "make" a nice design file once, instead of doing that awk magic each time. It seems wholly unnecessary
- Call the first file metadata.csv then out of that create design.csv and then the lines won't have excessive codes:
- Instead of:
```
make genome GENOME_ACC={}
```
put the actual accession number in there to show it works
- Then use the parallel construct:
```
cat design.csv | parallel --colsep , --header : make align SRR={SRR} SAMPLE={SAMPLE}
```
- Use --dry-run to see what parallel will do
---

As a reminder, my paper is titled: "Zika Virus Targets Human Cortical Neural Precursors and Attenuates Their Growth". This paper can be found at this link: https://pubmed.ncbi.nlm.nih.gov/26952870/. The acquistion number are as follows:
- BioProject: PRJNA313294
- GEO series: GSE78711
- SRR3194431

## Create a design.csv file that connects the SRR numbers to the sample names.
First, I got the metadata.
```
bio search PRJNA313294 -H --csv > metadata.csv
```

Then, I made my design file from the metadata. 
```
awk -F, 'BEGIN{OFS=","}
NR==1{
  for(i=1;i<=NF;i++) h[$i]=i
  print "Run,Sample,Layout,Description"; next
}
{
  run=$h["run_accession"]
  sample=$h["sample_alias"]
  layout=toupper($h["library_layout"])
  if(layout!="PE" && layout!="SE" && layout!="PAIRED" && layout!="SINGLE") layout="PE"
  if(layout=="PAIRED") layout="PE"
  if(layout=="SINGLE") layout="SE"
  description=$h["sample_description"]
  if(run!="" && sample!="") print run,sample,layout,description
}' metadata.csv | LC_ALL=C sort -u > design.csv
```
The result should look something like this:
```
Run,Sample,Layout,Description
SRR3191542,GSM2073121,PE,Mock1-1
SRR3191543,GSM2073122,PE,Mock2-1
SRR3191544,GSM2073123,PE,ZIKV1-1
SRR3191545,GSM2073124,PE,ZIKV2-1
SRR3194428,GSM2075585,SE,Mock1-2
SRR3194429,GSM2075586,SE,Mock2-2
SRR3194430,GSM2075587,SE,ZIKV1-2
SRR3194431,GSM2075588,SE,ZIKV2-2
```

## Usage
This pipeline processes RNA-Seq data for the Zika virus using a modular Makefile. Each step is defined as a separate target, allowing you to run them individually or all at once. The workflow supports both single-end and paired-end reads.

| Target        | Description                                                                 |
|---------------|-----------------------------------------------------------------------------|
| genome        | Downloads and extracts the Zika virus genome using NCBI Datasets CLI        |
| reads_single  | Downloads single-end reads from SRA and runs FastQC                         |
| reads_paired  | Downloads paired-end reads from SRA and runs FastQC                         |
| index         | Indexes the reference genome using bwa index                                |
| align_single  | Aligns single-end reads, sorts and indexes BAM, and generates alignment stats |
| align_paired  | Aligns paired-end reads, sorts and indexes BAM, and generates alignment stats |
| wiggle        | Converts BAM file to BigWig format using bedtools and bedGraphToBigWig      |
| clean         | Removes all generated files and folders (refs, reads, bam)                  |
| all           | Runs the full pipeline: genome, reads, index, align, wiggle                 |

Step-by-Step Usage
1. Download the Zika virus genome
```
make genome GENOME_ACC=GCF_000882815.3
```

2. Download sequencing reads
For single-end reads:
```
make reads_single SRR_ACC=SRR1972739 
```
For paired-end reads:
```
make reads_paired SRR_ACC=SRR1972740
```

3. Index the reference genome
```
make index REF=refs/zika_paper_genome.fa
```

4. Align reads to the genome
For single-end reads:
```
make align_single SRR_ACC=SRR1972739 REF=ref/zika_paper_genome.fa
```
For paired-end reads:
```
make align_paired SRR_ACC=SRR1972740 REF=ref/zika_paper_genome.fa
```

5. Convert BAM to BigWig
```
make wiggle SRR_ACC=SRR1972739 REF=ref/zika_paper_genome.fa
```

Extra: Clean Up
```
make clean
```

### Input & Output Summary
Input:
- Genome accession (GCF_000882815.3 for example)
- SRA accession (SRR1972739 for example)
- Reference genome: refs/zika_paper_genome.fa
Output:
- fastqc report
- Sorted BAM: bam/SRR1972739.bam for example
- Alignment report: bam/SRR1972739_alignment_report.txt for example
- BigWig: bam/SRR1972739.bw

## Running in parallel
### 1. Do the one time setup:
```
make genome GENOME_ACC=GCF_000882815.3
make index REF=refs/zika_paper_genome.fa
```

### 2. Parallelization 
First, using dry-run:
```
cat design.csv | parallel --dry-run --colsep , --header : '
  if [ "{Layout}" = "PE" ]; then
      make reads_paired align_paired wiggle SRA_ACC={Run} REF=refs/zika_paper_genome.fa;
  else
      make reads_single align_single wiggle SRA_ACC={Run} REF=refs/zika_paper_genome.fa;
  fi'
```

Now, onto the actual command:
```
cat design.csv | parallel -j 4 --colsep , --header : '
  if [ "{Layout}" = "PE" ]; then
      make reads_paired align_paired wiggle SRA_ACC={Run} REF=refs/zika_paper_genome.fa;
  else
      make reads_single align_single wiggle SRA_ACC={Run} REF=refs/zika_paper_genome.fa;
  fi'
```

