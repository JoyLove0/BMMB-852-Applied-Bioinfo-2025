# Week 11: Establish the effects of variants  
--
### Goals:
- Generate a VCF file.
- Evaluate the effects of at least three variants in your VCF file.
- You may use visual inspection via IGV, or you may use a software tool like snpEff, VEP, or any other alternative.
- Try to find variants with different types of effects.
- Summarize the process you followed and the results you obtained.
--

## Generate a VCF file
I used my Makefile to generate VCF files. To do this, I do the following:

```
make genome
make index
```
```
parallel -j4 --colsep ',' --header : 'if [ "{Layout}" = "PE" ]; then make reads_paired align_paired variants SRA_ACC={Run} REF=refs/zika_paper_genome.fa; else make reads_single align_single variants SRA_ACC={Run} REF=refs/zika_paper_genome.fa; fi' :::: design.csv
```

## Varaint Analysis
I chose the following file to do my analysis:
```SRR3194430_variants.vcf.gz```

Since the variants are sparse, I did the analysis visually using IGV.



The only type of variants are SNPs and they are all synonomous synonymous mutations. See below as an example: