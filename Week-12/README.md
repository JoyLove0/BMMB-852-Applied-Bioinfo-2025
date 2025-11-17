# Week 12: Evaluate data from the Cancer Genome in a Bottle project

I chose Option 2: Produce and evaluate variant calls for assignment. See details below:

Work with data from the Cancer Genome in a Bottle project and produce an evaluation. Generate variant calls from the Cancer Genome in a Bottle data and evaluate their quality.

# Tasks:

- Call variants for normal and tumor samples in a region of interest
- Compare variant calls between samples and identify tumor-specific variants
- Compare your results to the gold standard DeepVariant calls (if available)

## Downloading the Toolbox
```
bio code
```
Checking Installation:
```
cd src/run/
ls
```

Output:
```
aria.mk              bioproject.mk        curl.mk              freebayes.mk         ivar.mk              salmon.mk            star.mk              wiggle.mk
bcftools.mk          bowtie2.mk           datasets.mk          gatk.mk              minimap2.mk          snpeff.mk            tester.mk
bcftools_parallel.mk bwa.mk               deepvariant.mk       genbank.mk           refgenie.mk          splitchrom.mk        time.sh
bgzip.mk             coverage.mk          fastp.mk             hisat2.mk            rsync.mk             sra.mk               vep.mk
```

## Getting the Reference, Alignment, and VCF
```
bash get_ref_alignment_variants.sh
```

The code ```get_ref_alignment_variants.sh``` is as follows:
```
# ----------------------------------------------------------------------
# Set the URL of the reference genome
REF_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz

# Download the reference genome
REF=refs/GRCh38.fasta.gz

# The chromosome of interest.
CHR=chr12

# A subset of the reference genome for the region of interest.
REF=refs/${CHR}.fasta

# Create the reference directory
mkdir -p refs

# Download the reference genome
curl -L ${REF_URL} > ${REF}

# Index the reference genome for use with samtools.
samtools faidx ${REF}

# ----------------------------------------------------------------------
# Set the URL of the BAM file
BAM_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-N-D_Element-StdInsert_77x_GRCh38-GIABv3.bam

# The BAM file name.
BAM=bam/KRAS-N-P.bam

# Set the region of interest.
REGION=chr12:25,205,246-25,250,936

# Create the BAM directory
mkdir -p bam

# Extract the reads for the region of interest.
samtools view -b ${BAM_URL} ${REGION} > ${BAM}

# Index the BAM file.
samtools index ${BAM}

# Get statistics on the BAM file.
samtools flagstat ${BAM}
# ----------------------------------------------------------------------
VCF=vcf/KRAS-N-P.vcf.gz
# Run the variant calling for the region of interest.
make -f src/run/bcftools.mk \
        REF=${REF} \
        BAM=${BAM} \
        VCF=${VCF} \
        run

# Get the toolbox
bio code

VCF=vcf/KRAS-N-P.vcf.gz
# Run the variant calling for the region of interest.
make -f src/run/bcftools.mk \
        REF=${REF} \
        BAM=${BAM} \
        VCF=${VCF} \
        run

# ---------------------------------------------------------------------
```

Output:
```
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  754M  100  754M    0     0  36.5M      0  0:00:20  0:00:20 --:--:-- 19.1M
28035 + 0 in total (QC-passed reads + QC-failed reads)
27992 + 0 primary
0 + 0 secondary
43 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
27994 + 0 mapped (99.85% : N/A)
27951 + 0 primary mapped (99.85% : N/A)
27992 + 0 paired in sequencing
14001 + 0 read1
13991 + 0 read2
27783 + 0 properly paired (99.25% : N/A)
27910 + 0 with itself and mate mapped
41 + 0 singletons (0.15% : N/A)
113 + 0 with mate mapped to a different chr
113 + 0 with mate mapped to a different chr (mapQ>=5)
mkdir -p vcf/
bcftools mpileup -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP' -O u -f refs/chr12.fasta bam/KRAS-N-P.bam | \
        bcftools call --ploidy 2 --annotate 'FORMAT/GQ'  -mv -O u | \
        bcftools norm -f refs/chr12.fasta -d all -O u | \
        bcftools sort -O z > vcf/KRAS-N-P.vcf.gz
Writing to /var/folders/6g/hv_kp3790hl1j2kyx7_kf818sy3d__/T//bcftools.sneiWC
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 100
Note: The maximum per-sample depth with -d 100 is 100.0x
Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      79/0/0/17/0/0/0
Merging 1 temporary files
Done
Cleaning
bcftools index -t -f vcf/KRAS-N-P.vcf.gz
-rw-r--r--@ 1 jal7297  357253257   6.7K Nov 16 23:08 vcf/KRAS-N-P.vcf.gz
```

### Tumor
I reran this code for the tumor sample with the following changes:
```
BAM_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-T_Element-StdInsert_111x_GRCh38-GIABv3.bam
BAM=bam/KRAS-T.bam
VCF=vcf/KRAS-T.vcf.gz
```

Here is the new code:
```
# ----------------------------------------------------------------------
# Set the URL of the reference genome
REF_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz

# Download the reference genome
REF=refs/GRCh38.fasta.gz

# The chromosome of interest.
CHR=chr12

# A subset of the reference genome for the region of interest.
REF=refs/${CHR}.fasta

# Create the reference directory
mkdir -p refs

# Download the reference genome
curl -L ${REF_URL} > ${REF}

# Index the reference genome for use with samtools.
samtools faidx ${REF}

# ----------------------------------------------------------------------
# Set the URL of the BAM file
BAM_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-T_Element-StdInsert_111x_GRCh38-GIABv3.bam

# The BAM file name.
BAM=bam/KRAS-T.bam

# Set the region of interest.
REGION=chr12:25,205,246-25,250,936

# Create the BAM directory
mkdir -p bam

# Extract the reads for the region of interest.
samtools view -b ${BAM_URL} ${REGION} > ${BAM}

# Index the BAM file.
samtools index ${BAM}

# Get statistics on the BAM file.
samtools flagstat ${BAM}
# ----------------------------------------------------------------------
VCF=vcf/KRAS-N-P.vcf.gz
# Run the variant calling for the region of interest.
make -f src/run/bcftools.mk \
        REF=${REF} \
        BAM=${BAM} \
        VCF=${VCF} \
        run

# Get the toolbox
bio code

VCF=vcf/KRAS-T.vcf.gz
# Run the variant calling for the region of interest.
make -f src/run/bcftools.mk \
        REF=${REF} \
        BAM=${BAM} \
        VCF=${VCF} \
        run

# ----------------------------------------------------------------------
```

Ran:
```
bash get_ref_alignment_variants.sh 
```

Output:
```
 % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  754M  100  754M    0     0  66.5M      0  0:00:11  0:00:11 --:--:-- 68.1M
57327 + 0 in total (QC-passed reads + QC-failed reads)
57256 + 0 primary
0 + 0 secondary
71 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
57257 + 0 mapped (99.88% : N/A)
57186 + 0 primary mapped (99.88% : N/A)
57256 + 0 paired in sequencing
28626 + 0 read1
28630 + 0 read2
56816 + 0 properly paired (99.23% : N/A)
57116 + 0 with itself and mate mapped
70 + 0 singletons (0.12% : N/A)
282 + 0 with mate mapped to a different chr
279 + 0 with mate mapped to a different chr (mapQ>=5)
mkdir -p vcf/
bcftools mpileup -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP' -O u -f refs/chr12.fasta bam/KRAS-T.bam | \
        bcftools call --ploidy 2 --annotate 'FORMAT/GQ'  -mv -O u | \
        bcftools norm -f refs/chr12.fasta -d all -O u | \
        bcftools sort -O z > vcf/KRAS-N-P.vcf.gz
Writing to /var/folders/6g/hv_kp3790hl1j2kyx7_kf818sy3d__/T//bcftools.dFPpVG
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 100
Note: The maximum per-sample depth with -d 100 is 100.0x
Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      80/0/0/17/0/0/0
Merging 1 temporary files
Done
Cleaning
bcftools index -t -f vcf/KRAS-N-P.vcf.gz
-rw-r--r--@ 1 jal7297  357253257   6.3K Nov 16 23:17 vcf/KRAS-N-P.vcf.gz
# Skipped 69 files that already exist.
# Use -u to update existing files.
# Biostar Workflows: https://www.biostarhandbook.com/
mkdir -p vcf/
bcftools mpileup -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP' -O u -f refs/chr12.fasta bam/KRAS-T.bam | \
        bcftools call --ploidy 2 --annotate 'FORMAT/GQ'  -mv -O u | \
        bcftools norm -f refs/chr12.fasta -d all -O u | \
        bcftools sort -O z > vcf/KRAS-T.vcf.gz
[mpileup] 1 samples in 1 input files
Writing to /var/folders/6g/hv_kp3790hl1j2kyx7_kf818sy3d__/T//bcftools.oG18kh
[mpileup] maximum number of reads per input file set to -d 100
Note: The maximum per-sample depth with -d 100 is 100.0x
Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      80/0/0/17/0/0/0
Merging 1 temporary files
Done
Cleaning
bcftools index -t -f vcf/KRAS-T.vcf.gz
-rw-r--r--@ 1 jal7297  357253257   6.3K Nov 16 23:18 vcf/KRAS-T.vcf.gz
```

## Snpcalling

```
# Make the design file
make -f src/recipes/cancer-study-snpcall.mk design
# Batch process the samples in the design file
make -f src/recipes/cancer-study-snpcall.mk batch
```

Output:
```
# Created: design.csv
-rw-r--r--@ 1 jal7297  357253257   350B Nov 16 23:36 design.csv
cat design.csv | parallel \
        --header : --colsep , --lb \
        make -f src/recipes/cancer-study-snpcall.mk \
        NAME={name} \
        BAM_URL={url} \
        vcf
make[1]: Entering directory '/Users/jal7297/Desktop/Applied_Bioinfo/BMMB-852-Applied-Bioinfo-2025/Week-12'
make[1]: Entering directory '/Users/jal7297/Desktop/Applied_Bioinfo/BMMB-852-Applied-Bioinfo-2025/Week-12'
make -f src/run/bcftools.mk \
        REF=refs/GRCh38.fa.gz \
        BAM=bam/KRAS-T.bam \
        VCF=vcf/KRAS-T.vcf.gz \
        run
-rw-r--r--@ 1 jal7297  357253257   6.7K Nov 16 23:33 vcf/KRAS-N-P.vcf.gz
make[1]: Leaving directory '/Users/jal7297/Desktop/Applied_Bioinfo/BMMB-852-Applied-Bioinfo-2025/Week-12'
mkdir -p vcf/
bcftools mpileup -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP' -O u -f refs/GRCh38.fa.gz bam/KRAS-T.bam | \
        bcftools call --ploidy 2 --annotate 'FORMAT/GQ'  -mv -O u | \
        bcftools norm -f refs/GRCh38.fa.gz -d all -O u | \
        bcftools sort -O z > vcf/KRAS-T.vcf.gz
Writing to /var/folders/6g/hv_kp3790hl1j2kyx7_kf818sy3d__/T//bcftools.n28Yeg
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 100
Note: The maximum per-sample depth with -d 100 is 100.0x
bcftools index -t -f vcf/KRAS-T.vcf.gz
Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      80/0/0/17/0/0/0
Merging 1 temporary files
Done
Cleaning
-rw-r--r--@ 1 jal7297  357253257   6.3K Nov 16 23:36 vcf/KRAS-T.vcf.gz
-rw-r--r--@ 1 jal7297  357253257   6.3K Nov 16 23:36 vcf/KRAS-T.vcf.gz
make[1]: Leaving directory '/Users/jal7297/Desktop/Applied_Bioinfo/BMMB-852-Applied-Bioinfo-2025/Week-12'
```

## Variant Evaluation
Merge the various VCF files into a single VCF file and then find unqiue and common variants between the control and tumor samples.

### Merge VCF files:
```
bcftools merge -W -O z -o vcf/merged.vcf.gz vcf/KRAS-N-P.vcf.gz vcf/KRAS-T.vcf.gz
```

### Get the intersected variant calls:
```
bcftools isec -p vcf/isec vcf/KRAS-N-P.vcf.gz vcf/KRAS-T.vcf.gz
```

The output is as follows:
- vcf/isec/0000.vcf - variants unique to the control sample
- vcf/isec/0001.vcf - variants unique to the tumor sample
- vcf/isec/0002.vcf - variants present in both the control and tumor samples

### Renamed VCF files:
```
bcftools view -W -O z -o vcf/tumor.vcf.gz vcf/isec/0001.vcf
```

```
bcftools view -W -O z -o vcf/control.vcf.gz vcf/isec/0000.vcf
```

```
bcftools view -W -O z -o vcf/both.vcf.gz vcf/isec/0002.vcf
```

## Removing False Postitives:
```
bash false_positive.sh
```

Here is the code ```false_positive.sh```:
For Tumor:
```
# The statistics for the gold standard.
VCF=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz
bcftools view $VCF chr12:24000000-26000000 | bcftools stats

# The statistics for the currentvariant calls.
bcftools view vcf/tumor.vcf.gz | bcftools stats
```

Output:
```
# This file was produced by bcftools stats (1.22+htslib-1.22.1) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats 
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       <STDIN>
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      1
SN      0       number of records:      24
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 22
SN      0       number of MNPs: 0
SN      0       number of indels:       2
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
# TSTV, transitions/transversions
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)
#   - transversions, see https://en.wikipedia.org/wiki/Transversion
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       11      11      1.00    11      11      1.00
# SiS, Singleton stats:
#   - allele count, i.e. the number of singleton genotypes (AC=1)
#   - number of transitions, see above
#   - number of transversions, see above
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
SiS     0       1       22      11      11      2       0       0       2
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent    [10]not applicable
AF      0       0.000000        22      11      11      2       0       0       2
# QUAL, Stats by quality
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
QUAL    0       21.9    1       1       0       0
QUAL    0       22.6    2       0       2       0
QUAL    0       22.8    0       0       0       1
QUAL    0       24.5    0       0       0       1
QUAL    0       30.7    1       1       0       0
QUAL    0       33.5    1       0       1       0
QUAL    0       33.6    1       0       1       0
QUAL    0       33.7    1       1       0       0
QUAL    0       34.2    1       1       0       0
QUAL    0       34.6    1       0       1       0
QUAL    0       99.0    13      7       6       0
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]number of sites      [5]number of genotypes  [6]mean VAF
IDD     0       -20     1       0       .
IDD     0       -1      1       0       .
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     0
ST      0       A>G     0
ST      0       A>T     2
ST      0       C>A     5
ST      0       C>G     0
ST      0       C>T     4
ST      0       G>A     7
ST      0       G>C     2
ST      0       G>T     1
ST      0       T>A     1
ST      0       T>C     0
ST      0       T>G     0
# DP, depth:
#   - set id, see above
#   - the depth bin, corresponds to the depth (unless --depth was given)
#   - number of genotypes with this depth (zero unless -s/-S was given)
#   - fraction of genotypes with this depth (zero unless -s/-S was given)
#   - number of sites with this depth
#   - fraction of sites with this depth
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
# This file was produced by bcftools stats (1.22+htslib-1.22.1) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats 
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       <STDIN>
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      1
SN      0       number of records:      1
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
# TSTV, transitions/transversions
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)
#   - transversions, see https://en.wikipedia.org/wiki/Transversion
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       0       1       0.00    0       1       0.00
# SiS, Singleton stats:
#   - allele count, i.e. the number of singleton genotypes (AC=1)
#   - number of transitions, see above
#   - number of transversions, see above
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
SiS     0       1       1       0       1       0       0       0       0
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent    [10]not applicable
AF      0       0.000000        1       0       1       0       0       0       0
# QUAL, Stats by quality
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
QUAL    0       221.7   1       0       1       0
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]number of sites      [5]number of genotypes  [6]mean VAF
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     0
ST      0       A>G     0
ST      0       A>T     0
ST      0       C>A     1
ST      0       C>G     0
ST      0       C>T     0
ST      0       G>A     0
ST      0       G>C     0
ST      0       G>T     0
ST      0       T>A     0
ST      0       T>C     0
ST      0       T>G     0
# DP, depth:
#   - set id, see above
#   - the depth bin, corresponds to the depth (unless --depth was given)
#   - number of genotypes with this depth (zero unless -s/-S was given)
#   - fraction of genotypes with this depth (zero unless -s/-S was given)
#   - number of sites with this depth
#   - fraction of sites with this depth
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
DP      0       105     0       0.000000        1       100.000000
```

For control:
```
# The statistics for the gold standard.
VCF=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz
bcftools view $VCF chr12:24000000-26000000 | bcftools stats

# The statistics for the currentvariant calls.
bcftools view vcf/control.vcf.gz | bcftools stats
```

Output:
```
# This file was produced by bcftools stats (1.22+htslib-1.22.1) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats 
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       <STDIN>
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      1
SN      0       number of records:      24
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 22
SN      0       number of MNPs: 0
SN      0       number of indels:       2
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
# TSTV, transitions/transversions
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)
#   - transversions, see https://en.wikipedia.org/wiki/Transversion
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       11      11      1.00    11      11      1.00
# SiS, Singleton stats:
#   - allele count, i.e. the number of singleton genotypes (AC=1)
#   - number of transitions, see above
#   - number of transversions, see above
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
SiS     0       1       22      11      11      2       0       0       2
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent    [10]not applicable
AF      0       0.000000        22      11      11      2       0       0       2
# QUAL, Stats by quality
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
QUAL    0       21.9    1       1       0       0
QUAL    0       22.6    2       0       2       0
QUAL    0       22.8    0       0       0       1
QUAL    0       24.5    0       0       0       1
QUAL    0       30.7    1       1       0       0
QUAL    0       33.5    1       0       1       0
QUAL    0       33.6    1       0       1       0
QUAL    0       33.7    1       1       0       0
QUAL    0       34.2    1       1       0       0
QUAL    0       34.6    1       0       1       0
QUAL    0       99.0    13      7       6       0
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]number of sites      [5]number of genotypes  [6]mean VAF
IDD     0       -20     1       0       .
IDD     0       -1      1       0       .
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     0
ST      0       A>G     0
ST      0       A>T     2
ST      0       C>A     5
ST      0       C>G     0
ST      0       C>T     4
ST      0       G>A     7
ST      0       G>C     2
ST      0       G>T     1
ST      0       T>A     1
ST      0       T>C     0
ST      0       T>G     0
# DP, depth:
#   - set id, see above
#   - the depth bin, corresponds to the depth (unless --depth was given)
#   - number of genotypes with this depth (zero unless -s/-S was given)
#   - fraction of genotypes with this depth (zero unless -s/-S was given)
#   - number of sites with this depth
#   - fraction of sites with this depth
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
# This file was produced by bcftools stats (1.22+htslib-1.22.1) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats 
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       <STDIN>
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      1
SN      0       number of records:      0
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 0
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
# TSTV, transitions/transversions
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)
#   - transversions, see https://en.wikipedia.org/wiki/Transversion
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       0       0       0.00    0       0       0.00
# SiS, Singleton stats:
#   - allele count, i.e. the number of singleton genotypes (AC=1)
#   - number of transitions, see above
#   - number of transversions, see above
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
SiS     0       1       0       0       0       0       0       0       0
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent    [10]not applicable
# QUAL, Stats by quality
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]number of sites      [5]number of genotypes  [6]mean VAF
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     0
ST      0       A>G     0
ST      0       A>T     0
ST      0       C>A     0
ST      0       C>G     0
ST      0       C>T     0
ST      0       G>A     0
ST      0       G>C     0
ST      0       G>T     0
ST      0       T>A     0
ST      0       T>C     0
ST      0       T>G     0
# DP, depth:
#   - set id, see above
#   - the depth bin, corresponds to the depth (unless --depth was given)
#   - number of genotypes with this depth (zero unless -s/-S was given)
#   - fraction of genotypes with this depth (zero unless -s/-S was given)
#   - number of sites with this depth
#   - fraction of sites with this depth
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
```

For both:
```
# The statistics for the gold standard.
VCF=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz
bcftools view $VCF chr12:24000000-26000000 | bcftools stats

# The statistics for the currentvariant calls.
bcftools view vcf/both.vcf.gz | bcftools stats
```

Output:
```
# This file was produced by bcftools stats (1.22+htslib-1.22.1) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats 
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       <STDIN>
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      1
SN      0       number of records:      24
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 22
SN      0       number of MNPs: 0
SN      0       number of indels:       2
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
# TSTV, transitions/transversions
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)
#   - transversions, see https://en.wikipedia.org/wiki/Transversion
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       11      11      1.00    11      11      1.00
# SiS, Singleton stats:
#   - allele count, i.e. the number of singleton genotypes (AC=1)
#   - number of transitions, see above
#   - number of transversions, see above
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
SiS     0       1       22      11      11      2       0       0       2
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent    [10]not applicable
AF      0       0.000000        22      11      11      2       0       0       2
# QUAL, Stats by quality
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
QUAL    0       21.9    1       1       0       0
QUAL    0       22.6    2       0       2       0
QUAL    0       22.8    0       0       0       1
QUAL    0       24.5    0       0       0       1
QUAL    0       30.7    1       1       0       0
QUAL    0       33.5    1       0       1       0
QUAL    0       33.6    1       0       1       0
QUAL    0       33.7    1       1       0       0
QUAL    0       34.2    1       1       0       0
QUAL    0       34.6    1       0       1       0
QUAL    0       99.0    13      7       6       0
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]number of sites      [5]number of genotypes  [6]mean VAF
IDD     0       -20     1       0       .
IDD     0       -1      1       0       .
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     0
ST      0       A>G     0
ST      0       A>T     2
ST      0       C>A     5
ST      0       C>G     0
ST      0       C>T     4
ST      0       G>A     7
ST      0       G>C     2
ST      0       G>T     1
ST      0       T>A     1
ST      0       T>C     0
ST      0       T>G     0
# DP, depth:
#   - set id, see above
#   - the depth bin, corresponds to the depth (unless --depth was given)
#   - number of genotypes with this depth (zero unless -s/-S was given)
#   - fraction of genotypes with this depth (zero unless -s/-S was given)
#   - number of sites with this depth
#   - fraction of sites with this depth
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
# This file was produced by bcftools stats (1.22+htslib-1.22.1) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats 
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       <STDIN>
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      1
SN      0       number of records:      79
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 60
SN      0       number of MNPs: 0
SN      0       number of indels:       19
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
# TSTV, transitions/transversions
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)
#   - transversions, see https://en.wikipedia.org/wiki/Transversion
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       40      20      2.00    40      20      2.00
# SiS, Singleton stats:
#   - allele count, i.e. the number of singleton genotypes (AC=1)
#   - number of transitions, see above
#   - number of transversions, see above
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
SiS     0       1       7       3       4       7       0       0       7
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent    [10]not applicable
AF      0       0.000000        7       3       4       7       0       0       7
AF      0       0.990000        53      37      16      12      0       0       12
# QUAL, Stats by quality
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
QUAL    0       82.3    0       0       0       1
QUAL    0       122.8   0       0       0       1
QUAL    0       186.9   0       0       0       1
QUAL    0       210.1   0       0       0       1
QUAL    0       218.2   0       0       0       1
QUAL    0       218.4   1       0       1       0
QUAL    0       219.0   0       0       0       1
QUAL    0       219.9   0       0       0       1
QUAL    0       221.0   0       0       0       1
QUAL    0       222.1   0       0       0       1
QUAL    0       222.3   2       2       0       0
QUAL    0       222.4   4       1       3       0
QUAL    0       225.4   36      25      11      0
QUAL    0       228.0   0       0       0       2
QUAL    0       228.1   0       0       0       3
QUAL    0       228.2   0       0       0       3
QUAL    0       228.3   17      12      5       1
QUAL    0       228.4   0       0       0       1
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]number of sites      [5]number of genotypes  [6]mean VAF
IDD     0       -26     1       0       .
IDD     0       -5      1       0       .
IDD     0       -4      2       0       .
IDD     0       -3      1       0       .
IDD     0       -2      2       0       .
IDD     0       -1      5       0       .
IDD     0       1       5       0       .
IDD     0       5       1       0       .
IDD     0       13      1       0       .
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     4
ST      0       A>G     5
ST      0       A>T     4
ST      0       C>A     3
ST      0       C>G     2
ST      0       C>T     14
ST      0       G>A     13
ST      0       G>C     1
ST      0       G>T     1
ST      0       T>A     2
ST      0       T>C     8
ST      0       T>G     3
# DP, depth:
#   - set id, see above
#   - the depth bin, corresponds to the depth (unless --depth was given)
#   - number of genotypes with this depth (zero unless -s/-S was given)
#   - fraction of genotypes with this depth (zero unless -s/-S was given)
#   - number of sites with this depth
#   - fraction of sites with this depth
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
DP      0       46      0       0.000000        1       1.265823
DP      0       47      0       0.000000        1       1.265823
DP      0       52      0       0.000000        1       1.265823
DP      0       53      0       0.000000        2       2.531646
DP      0       54      0       0.000000        1       1.265823
DP      0       61      0       0.000000        1       1.265823
DP      0       63      0       0.000000        2       2.531646
DP      0       64      0       0.000000        2       2.531646
DP      0       65      0       0.000000        2       2.531646
DP      0       66      0       0.000000        1       1.265823
DP      0       68      0       0.000000        1       1.265823
DP      0       69      0       0.000000        3       3.797468
DP      0       71      0       0.000000        1       1.265823
DP      0       73      0       0.000000        2       2.531646
DP      0       74      0       0.000000        1       1.265823
DP      0       75      0       0.000000        2       2.531646
DP      0       76      0       0.000000        1       1.265823
DP      0       77      0       0.000000        1       1.265823
DP      0       78      0       0.000000        2       2.531646
DP      0       79      0       0.000000        3       3.797468
DP      0       80      0       0.000000        2       2.531646
DP      0       81      0       0.000000        2       2.531646
DP      0       83      0       0.000000        1       1.265823
DP      0       84      0       0.000000        5       6.329114
DP      0       85      0       0.000000        1       1.265823
DP      0       86      0       0.000000        3       3.797468
DP      0       87      0       0.000000        2       2.531646
DP      0       88      0       0.000000        5       6.329114
DP      0       89      0       0.000000        2       2.531646
DP      0       90      0       0.000000        2       2.531646
DP      0       91      0       0.000000        3       3.797468
DP      0       92      0       0.000000        1       1.265823
DP      0       93      0       0.000000        1       1.265823
DP      0       95      0       0.000000        1       1.265823
DP      0       96      0       0.000000        2       2.531646
DP      0       97      0       0.000000        2       2.531646
DP      0       98      0       0.000000        2       2.531646
DP      0       99      0       0.000000        2       2.531646
DP      0       100     0       0.000000        3       3.797468
DP      0       101     0       0.000000        4       5.063291
DP      0       102     0       0.000000        2       2.531646
```

### Filtering Variants

Tumor:
```
bcftools view --max-alleles 2 --types snps -O z -o vcf/tumor.filtered.vcf.gz vcf/tumor.vcf.gz
```

Control:
```
bcftools view --max-alleles 2 --types snps -O z -o vcf/control.filtered.vcf.gz vcf/control.vcf.gz
```

Both:
```
bcftools view --max-alleles 2 --types snps -O z -o vcf/both.filtered.vcf.gz vcf/both.vcf.gz
```

For note, the Gold Standard is below:
```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz
```

# Deliverable: 
Report with variant counts, tumor-specific variants, and comparison to gold standard.

## Tumor Counts
### Gold Standard Variant Set (chr12:24–26 Mb)
| Metric        | Count  |
| ------------- | ------ |
| Total records | **24** |
| SNPs          | **22** |
| Indels        | **2**  |


### Tumor Variant Counts
| Metric        | Count |
| ------------- | ----- |
| Total records | **1** |
| SNPs          | **1** |
| Indels        | **0** |

### Control Variant Counts
| Metric        | Count |
| ------------- | ----- |
| Total records | **0** |
| SNPs          | **0** |
| Indels        | **0** |

### Summary Table
### Variant Calling Summary

| Dataset        | SNPs | Indels | Total | Key Notes |
|----------------|------|--------|--------|-----------|
| Gold standard  | 22   | 2      | 24     | Reference truth set |
| Tumor          | 1    | 0      | 1      | Severe under-calling; missing most true variants |
| Control        | 0    | 0      | 0      | No variants detected |
| Both samples   | 60   | 19     | 79     | Detects all true variants but adds ~55 false positives |

In comparison to the Golden Standard, both the tumor and control variant counts. The tumor sample contains only a single SNP, indicating substantial under-calling relative to the expected 24 true variants. The control sample contains no detected variants. The combined tumor–control callset detects all gold-standard variants but introduces a large number of false positives.
