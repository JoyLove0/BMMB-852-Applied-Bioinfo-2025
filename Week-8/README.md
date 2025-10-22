# Week 8: Automate a large scale analysis

## This assignment requires the presence of a Makefile, a README.md markdown file, and a design.csv file.
All files can be found in the Week-8 repo.

## 1. Identify the sample names that connect the SRR numbers to the samples.
As a reminder, my paper is titled: "Zika Virus Targets Human Cortical Neural Precursors and Attenuates Their Growth". This paper can be found at this link: https://pubmed.ncbi.nlm.nih.gov/26952870/. The acquistion number is as follows:
- BioProject: PRJNA313294
- GEO series: GSE78711
- SRR3194431

## 2. Create a design.csv file that connects the SRR numbers to the sample names.
Usually to do this, I would do as follows:
```
bio search PRJNA313294 -H --csv > design.csv
```

However, I got this error:
```
Error for https://www.ebi.ac.uk/ena/portal/api/filereport, {'accession': 'PRJNA31329', 'fields': 'run_accession,sample_accession,sample_alias,sample_description,first_public,country,scientific_name,fastq_bytes,base_count,read_count,library_name,library_strategy,library_source,library_layout,instrument_platform,instrument_model,study_title,fastq_ftp', 'result': 'read_run'}: 500 Server Error: Internal Server Error for url: https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA31329&fields=run_accession%2Csample_accession%2Csample_alias%2Csample_description%2Cfirst_public%2Ccountry%2Cscientific_name%2Cfastq_bytes%2Cbase_count%2Cread_count%2Clibrary_name%2Clibrary_strategy%2Clibrary_source%2Clibrary_layout%2Cinstrument_platform%2Cinstrument_model%2Cstudy_title%2Cfastq_ftp&result=read_run
```
This error suggests that the website that was storing this project no longer exsists. When I went to the url provided in the error message, I got a 404 error.
```
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡀⣤⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⠀⢀⣄⡀⠀⠀⢀⣴⣾⣿⣿⣿⣿⣿⣷⣦⣄⡀⠀⠀⣤⣟⠛⠋⠙⢷⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡴⠿⣫⠟⠛⢋⣧⣤⣶⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣼⡯⣉⠉⠓⢦⣀⠉⢢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⡴⠋⠀⣠⠖⠉⢋⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣄⠑⢄⠀⠈⢳⡀⣧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⢿⠏⠀⡴⠊⠁⢀⣴⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⢿⢦⣀⢳⠀⠀⢣⠸⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⡞⢀⠞⠀⠀⣴⣿⡿⡿⢿⣿⡟⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⣄⣀⠀⠉⠃⠀⠀⠳⢽⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣾⠏⣠⠋⠀⢀⡴⠛⠁⢀⣨⣾⡿⣇⠀⠛⠛⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢲⡀⠀⠀⠀⠀⠹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠏⠁⠀⠀⠀⠴⠋⠀⢀⣴⣿⣿⠟⠀⠈⣰⣦⣀⢀⣈⣿⠟⡿⢿⣿⠏⣈⣟⢻⣿⣿⣿⡀⠳⡄⠀⠀⠀⠀⠙⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠏⠀⠀⠀⠀⠀⠀⢀⣴⣿⣿⣿⠃⠀⢀⡄⠉⣉⠿⠛⠧⣞⣀⣧⡜⢛⣛⣻⣏⣀⢻⣿⣿⠃⠣⡇⠀⠀⠀⠀⠀⠘⢆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠏⠀⠀⠀⠀⠀⠀⢠⣾⣿⣿⣿⡏⠀⠀⠎⠀⠽⢟⣻⣟⡆⣸⠆⠸⡄⢻⣛⣛⠃⠘⢸⣿⡿⠀⠀⠸⡀⠀⠀⠀⠀⠀⠈⠳⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⣠⠏⠀⠀⠀⠀⠀⠀⢀⡞⢹⣿⣿⣿⡇⠀⠀⠀⠀⠀⠈⢁⣠⢤⡁⠀⠀⢈⢦⡀⠀⠀⠀⢸⣿⠇⠀⠀⠀⠳⣄⠀⠀⠀⠀⠀⠀⠙⢦⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⣠⠊⠁⠀⠀⠀⠀⠀⠀⢀⡞⠀⠈⣿⠙⠛⣇⠀⠀⠀⠀⠀⠀⠉⢀⡼⠇⠀⠀⠘⢦⡁⠀⠀⠀⢸⣟⠀⠀⠀⠀⠀⠀⠙⢦⠀⠀⠀⠀⠀⠀⠱⣄⠀⠀⠀⠀⠀⠀
⠀⠀⠀⡴⣯⡟⠁⠀⠀⠀⠀⠀⢀⡤⠖⠋⠀⠀⠀⠸⡎⠿⠏⠀⠀⠀⠀⠀⠀⠀⠚⠓⠚⠣⠴⠚⠛⣃⠀⠀⠀⣺⡏⠀⠀⠀⠀⠀⠀⠀⠈⢣⠀⠀⠀⠀⠀⠀⠈⠳⠀⠀⠀⠀⠀
⠀⠀⠘⠀⢈⡏⠀⠀⠀⠀⠀⢀⡞⠀⠀⠀⠀⠀⠀⠀⢻⣄⣠⡄⠀⠀⠀⠀⠀⠀⣴⠋⣠⣤⣤⣄⠀⠘⠆⠀⠀⡟⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⢣⣀⣀⣀⣀⣤⣤⣤⠬⠖⠀⠀⠀
⠤⠀⠀⠀⢸⡀⠀⠀⠀⠀⢠⠎⠀⠀⠀⠀⠀⠀⠀⠀⠘⣿⢿⡇⠀⠀⠀⠀⠀⠀⠁⠾⣗⠓⠒⢚⡗⠀⠀⠀⡼⠙⣆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠱
⠀⠀⠀⠀⠀⠉⠓⠒⠲⢰⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣤⡎⢣⡀⠀⠀⠀⠀⠀⠀⠀⠀⢉⣉⠁⠀⠀⠀⣼⠁⠀⠈⠳⣄⢀⡴⠒⠦⣄⣀⣀⣀⣽⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⣧⣄⣀⡤⠶⠤⣴⠒⠒⠒⠋⠁⢸⠁⠀⠙⢦⡀⠀⠀⠀⠀⠀⠘⠉⠈⠙⠂⢀⠞⢸⡆⠀⠀⠀⡼⠉⠀⠀⠀⠀⠀⠁⠀⠘⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⠀⠀⠀⠀⠀⠈⠙⠂⠀⠀⠀⢻⡀⠀⠀⠀⠈⠻⡒⠦⠤⢄⣀⣀⠀⣀⠴⠋⢀⡼⢧⠀⠀⣸⠁⠀⠀⠀⠀⠸⡀⠀⠀⠀⢳⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⡿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⢧⡀⠀⠀⠀⠀⣉⣀⣀⠴⢫⠇⠀⠙⠲⠧⡄⠀⠀⠀⠀⠀⢣⠀⠀⠀⠈⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⢀⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠲⢤⣀⣸⡉⠉⠉⠉⠉⠉⠀⠀⠀⡼⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠘⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠇⠀⠀⠀⠀⠀⠀⠀⢠⡇⠀⠀⠀⠀⢰⠇⠀⠀⠀⠀⠀⠀⢹⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
```
However, the show must go on! So, I went ahead and just downloaded the metadata for this project using the SRA Run Selector hosted by the NCBI. Here is the link: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA313294&o=acc_s%3Aa. With that downloaded file, I was able to make my design file with:
```
mv SraRunTable.csv design.csv
```

## 3. Create a Makefile that can produce multiple BAM alignment files (you can reuse the one from the previous assignment) where from an SRR you can produce a BAM alignment file named after the sample name.
This files is in the Week-8 repo.

## 4. Using GNU parallel run the Makefile on all (or at least 10) samples.
I had to run this twice. Once for the single reads and again for the paired-end reads. Also, in order to not run the genome acquistion step multiple times, the ``` make genome ``` needs to be ran seperate from the rest. 
```
make genome GENOME_ACC=GCF_000882815.3
```
```
awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="LibraryLayout") col=i} NR>1 && $col=="SINGLE" {print $1}' design.csv | \
parallel --jobs 4 --bar '
  make index REF=refs/zika_paper_genome.fa GENOME_ACC={BioProject};
  make reads_single SRA_ACC={Run};
  make align_single SRA_ACC={Run} REF=refs/zika_paper_genome.fa GENOME_ACC={BioProject};
  make wiggle SRA_ACC={} REF=refs/zika_paper_genome.fa
```
```
awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="LibraryLayout") col=i} NR>1 && $col=="PAIRED" {print $1}' design.csv | \
parallel --jobs 4 --bar '
  make index REF=refs/zika_paper_genome.fa GENOME_ACC={};
  make reads_paired SRA_ACC={};
  make align_paired SRA_ACC={} REF=refs/zika_paper_genome.fa GENOME_ACC={};
  make wiggle SRA_ACC={} REF=refs/zika_paper_genome.fa
```
## 5. Create a README.md file that explains how to run the Makefile/ Your README.md file should explain how to run the Makefile on the design.csv file.
### The result should consist of:
### - FASTQ read data named by the samples
### - FASTQC reports for each read
### - Alignments and coverage files in BAM and BW formats.
### - A statistics alignment report for BAM file

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
make genome GENOME_ACC={}
```

2. Download sequencing reads
For single-end reads:
```
make reads_single SRA_ACC={}
```
For paired-end reads:
```
make reads_paired SRA_ACC={}
```

3. Index the reference genome
```
make index REF=refs/zika_paper_genome.fa
```

4. Align reads to the genome
For single-end reads:
```
make align_single SRA_ACC={} REF=refs/zika_paper_genome.fa GENOME_ACC={}

```
For paired-end reads:
```
make align_paired SRA_ACC={} REF=refs/zika_paper_genome.fa GENOME_ACC={}
```

5. Convert BAM to BigWig
```
make wiggle SRA_ACC={} REF=refs/zika_paper_genome.fa

```

Extra: Clean Up
```
make clean
```

## Input & Output Summary
Input:
- Genome accession (GCF_000882815.3 for example)
- SRA accession (SRR1972739 for example)
- Reference genome: refs/zika_paper_genome.fa
Output:
- fastqc report
- Sorted BAM: bam/SRR1972739.bam for example
- Alignment report: bam/SRR1972739_alignment_report.txt for example
- BigWig: bam/SRR1972739.bw
