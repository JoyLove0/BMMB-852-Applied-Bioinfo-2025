# Week 10: Generate a multisample variant call file (VCF)
---
# Added Features
### Extend your existing Makefile to call variants on a single sample of your choice.
- Use your existing BAM file as input
- Generate a VCF file for this sample
- Follow best practices for variant calling
- Visualize the resulting VCF file data alongside the BAM file

### Call variants for all samples
- Run the variant calling workflow for all samples using your design.csv file.

### Create a multisample VCF
- Merge all individual sample VCF files into a single multisample VCF file (hint: bcftools merge)
- Visualize the multisample VCF in the context of the GFF annotation file.
---

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

| Target              | Description                                                                                                                                                                                        |
| :------------------ | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **genome**          | Downloads and extracts the Zika virus genome using the NCBI Datasets CLI                                                                                                                           |
| **reads_single**    | Downloads single-end reads from SRA and runs FastQC quality control                                                                                                                                |
| **reads_paired**    | Downloads paired-end reads from SRA and runs FastQC quality control                                                                                                                                |
| **index**           | Indexes the reference genome using `bwa index`                                                                                                                                                     |
| **align_single**    | Aligns single-end reads, sorts and indexes BAM, and generates alignment statistics                                                                                                                 |
| **align_paired**    | Aligns paired-end reads, sorts and indexes BAM, and generates alignment statistics                                                                                                                 |
| **wiggle**          | Converts BAM files to BigWig format using `bedtools genomecov` and `bedGraphToBigWig`                                                                                                              |
| **variants_sample** | Calls variants for a single sample using `bcftools mpileup` and `bcftools call`, indexes the resulting VCF, and automatically loads genome, BAM, GFF, and VCF tracks into IGV via its network port |
| **variants_multi**  | Merges multiple individual sample VCFs into a multisample VCF using `bcftools merge`, indexes it, and automatically loads all BAMs and the multisample VCF into IGV                                |
| **clean**           | Removes all generated files and folders (`refs`, `reads`, `bam`)                                                                                                                                   |
| **all**             | Runs the full pipeline: `genome`, `reads`, `index`, `align`, and `wiggle`                                                                                                                          |

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

6. Call Variants for a Single Sample
This step performs variant calling following bcftools best practices, producing a compressed and indexed .vcf.gz file for one sample.
```
make variants SRR_ACC=SRR1972739 REF=ref/zika_paper_genome.fa  F1='-d 100 --annotate ...' F2='--ploidy 2 --annotate ...' IGV_PORT=60151
```

7. Merge VCFs and Visualize Multisample Data
This step merges multiple sample VCFs into one multisample VCF, indexes it, and automatically visualizes all data in IGV.
```
ake variants_multi IGV_PORT=60151 VCF_LIST='vcf/SRR1.vcf.gz vcf/SRR2.vcf.gz'
```

Extra: Clean Up
```
make clean
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

