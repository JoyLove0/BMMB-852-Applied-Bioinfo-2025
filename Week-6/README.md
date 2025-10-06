# Week 6: Generate a BAM alignment file

## Make Makefile that includes rules for:
- Obtaining the genome
- Downloading sequencing reads from SRA
- index: Index the genome
- align: Generate a sorted and indexed BAM file by aligning reads to the genome

This file is in the Week-6 directory. 

## Explain the use of the Makefile in your project.

### Overview of Zika Virus Genome Alignment Pipeline

This project automates a basic bioinformatics workflow using a **Makefile**. It downloads the Zika virus genome and sequencing reads, indexes the genome, and aligns the reads to produce a sorted and indexed BAM file.

The pipeline includes the following steps:

1. **Obtain the reference genome** (Zika virus)
2. **Download sequencing reads** from the SRA
3. **Index the genome** using `bwa`
4. **Align reads** to the genome, sort, and index the output BAM

Each of these steps is defined as a separate **Makefile target**, allowing you to run them individually or all at once.

| **Target** | **Description**                                          |
| ---------- | -------------------------------------------------------- |
| `genome`   | Downloads and extracts the Zika virus genome             |
| `reads`    | Downloads reads from SRA                                 |
| `index`    | Indexes the reference genome using `bwa index`           |
| `align`    | Aligns reads, sorts SAM to BAM, and indexes the BAM file |
| `clean`    | Removes all generated files and folders                  |
| `all`      | Default target, runs everything up to `align`            |

### Usage
You can run the entire pipeline or individual steps using the make command.

### Run the Entire Pipeline
```
make
```
This is equivalent to:
```
make align
```
It will:
- Download the genome
- Download reads
- Index the genome
- Align reads and produce sorted + indexed BAM

### Step-by-Step Usage
1. ``` make genomen ```
Downloads and extracts the Zika virus genome.
Source: NCBI Genome via datasets CLI
Output: refs/zika_paper_genome.fa

2.``` make reads ```
Downloads sequencing reads from the SRA.
Accession: SRR9945583
Output: reads/SRR9945583_1.fastq
Only single-end reads are used in this example.

3. ```make index```
Indexes the reference genome using BWA.
Input: refs/zika_paper_genome.fa
Output: BWA index files (.bwt, .pac, etc.)
Note: make index depends on make genome.

4. ```make align```
Aligns reads to the genome, sorts the result, and creates an indexed BAM.
Input: Genome, reads
Output: bam/SRR9945583.bam and bam/SRR9945583.bam.bai
Note: make align depends on make index and make reads.

Extra: Clean Up 
``` make clean ```
Deletes all generated directories and files:
- refs/
- reads/
- bam/
  
### Customization
You can change the genome or read accession numbers by editing the top of the Makefile:
```
GENOME_ACC = GCF_000882815.3
SRA_ACC = SRR9945583
```
To process another sample, update SRA_ACC and re-run the relevant targets.

## Visualize the resulting BAM files for both simulated reads and reads downloaded from SRA.
My SRR file aligned to nothing at first (see below). 

<img width="1440" height="900" alt="Screenshot 2025-10-05 at 11 17 43 PM" src="https://github.com/user-attachments/assets/38d5e4c7-ba5a-40d5-99ba-610209eecf4d" />

So, I had to download the whole SRA file to see if I got any hits. I did! But the hit are very suspicious. Note that the IGV idenifies paired read but this is single read data. Also, see the unusually high covergae on the right end.

<img width="1440" height="900" alt="Screenshot 2025-10-05 at 11 51 59 PM" src="https://github.com/user-attachments/assets/69016acf-87f4-4b2e-8ce4-409a3ab3bbbe" />

The next section will emphasize these suspicions. 

## Generate alignment statistics for the BAM file.
- What percentage of reads aligned to the genome?

0% of the reads aligned to the reference genone, even after downloading the entire SRR file.  ¯\_(ツ)_/¯     

Code:
```
samtools flagstat bam/SRR9945583.bam
```

Output:
```
26381611 + 0 in total (QC-passed reads + QC-failed reads)
26381611 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
195 + 0 mapped (0.00% : N/A)
195 + 0 primary mapped (0.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)```
```
Also, when I copied the details from the covergae track on IGV, I got nosense:
```
NC_012532.1:10,717
<hr>Total count: 157
A      : 0
C      : 157  (100%,     157+,   0- )
G      : 0
T      : 0
N      : 0
```

- What was the expected average coverage?

For the ~1500 reads I used orginally, I should have got 10X coverage. When using all reads (10,794 bp), I should have got 464,041X coverage.

How 464,041X:

coverage = total bases/ genome size = 5,008,831,278/ 10,794 = 464,041.

- What is the observed average coverage?

0X  ｡°(°¯᷄◠¯᷅°)°｡

- How much does the coverage vary across the genome? (Provide a visual estimate.)

The coverage is practically non-existant expect for at a one point at the end and that is likely due to error. 

<img width="1429" height="147" alt="Screenshot 2025-10-06 at 12 11 43 AM" src="https://github.com/user-attachments/assets/83b3fe3c-91dc-4b64-b19f-9ea0298f6666" />

---
In conclusion:


![aZrXrgQ_700bwp](https://github.com/user-attachments/assets/d8b1eaa6-db09-4a75-b0a9-b6c5ea2d2c86)


