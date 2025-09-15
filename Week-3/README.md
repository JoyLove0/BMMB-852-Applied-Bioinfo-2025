# Assigment 3

For this assignment, the NCBI's reference genomes for the Black Rhino (Diceros bicornis) was used. The link to this can be found at: [Black Rhino Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_020826845.1/). This information can also be found [here](https://hgdownload.soe.ucsc.edu/hubs/mammals/index.html). Once at that sight navigate to the same refernce genome: GCF_020826845.1_mDicBic1.mat.cur (number 89). I am noting this because this information is more accessible in IGV. 

The genome and annotation follow were acquired as follows:

Command:
```
$ curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/826/845/GCF_020826845.1_mDicBic1.mat.cur/GCF_020826845.1_mDicBic1.mat.cur_genomic.fna.gz
$ gunzip GCF_020826845.1_mDicBic1.mat.cur_genomic.fna.gz 
```
Output:
```
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  858M  100  858M    0     0  5198k      0  0:02:49  0:02:49 --:--:-- 6580k
```

Command:
```
$ curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/826/845/GCF_020826845.1_mDicBic1.mat.cur/GCF_020826845.1_mDicBic1.mat.cur_genomic.gff.gz
$ gunzip GCF_020826845.1_mDicBic1.mat.cur_genomic.gff.gz
```
Output:
```
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 14.6M  100 14.6M    0     0  9735k      0  0:00:01  0:00:01 --:--:-- 9729k
```

## Use IGV to visualize your genome and the annotations relative to the genome.

I found the best way to view the genome and annotations files is through IGV itself. See below:

1. Navigated to Genome > Load Genome from UCSC GenArk...
2. Searched GCF_020826845.1 and selected it. Hit load
3. Used default selectionions


## How big is the genome, and how many features of each type does the GFF file contain?
Command:
```
$ awk '{sum += $2} END {print "Total genome size (bp):", sum}' GCF_020826845.1_mDicBic1.mat.cur_genomic.fna.fai

```
Output:
```
Total genome size (bp): 3005535620
```

Command:
```
$ grep -v '^#' GCF_020826845.1_mDicBic1.mat.cur_genomic.gff | cut -f3 | sort | uniq -c | sort -nr
```
Output:
```
648284 exon
569420 CDS
47787 mRNA
27017 gene
7612 cDNA_match
5465 lnc_RNA
3351 pseudogene
1600 transcript
1076 region
 579 tRNA
 537 snoRNA
 487 snRNA
 371 rRNA
 142 V_gene_segment
  24 C_gene_segment
  23 ncRNA
```

The genome is 3005535620 bp and the features with its corresponding counts is as follows: 

| Feature         | Count   |
|-----------------|---------|
| exon            | 648,284 |
| CDS             | 569,420 |
| mRNA            | 47,787  |
| gene            | 27,017  |
| cDNA_match      | 7,612   |
| lnc_RNA         | 5,465   |
| pseudogene      | 3,351   |
| transcript      | 1,600   |
| region          | 1,076   |
| tRNA            | 579     |
| snoRNA          | 537     |
| snRNA           | 487     |
| rRNA            | 371     |
| V_gene_segment  | 142     |
| C_gene_segment  | 24      |
| ncRNA           | 23      |

## From your GFF file, separate the intervals of type "gene" or "transcript" into a different file. Show the commands you used to do this.

Command:
```
$ awk '$3 == "gene" || $3 == "transcript"' GCF_020826845.1_mDicBic1.mat.cur_genomic.gff > black_rhino_genes_and_transcripts.gff
```
Output:
Snipett of file
```
NC_080740.1     Gnomon  gene    1420918 1438648 .       + $
NC_080740.1     Gnomon  gene    1448060 1483323 .       + $
NC_080740.1     Gnomon  gene    1488180 1508153 .       - $
NC_080740.1     Gnomon  gene    1520423 1608183 .       + $
NC_080740.1     Gnomon  gene    1615707 1677711 .       + $
NC_080740.1     Gnomon  transcript      1615707 1669947 . $
NC_080740.1     Gnomon  transcript      1615707 1669943 . $
NC_080740.1     Gnomon  transcript      1615707 1669943 . $

```
## Visualize the simplified GFF in IGV as a separate track. Compare the visualization of the original GFF with the simplified GFF.

The following command made the gene and transcrpt files into a track:
```
$ sort -k1,1 -k4,4n black_rhino_genes_and_transcripts.gff > black_rhino_genes_and_transcripts.sorted.gff
$ bgzip black_rhino_genes_and_transcripts.sorted.gff
$ tabix -p gff black_rhino_genes_and_transcripts.sorted.gff.gz
```

## Zoom in to see the sequences, expand the view to show the translation table in IGV. Note how the translation table needs to be displayed in the correct orientation for it to make sense.

## Visually verify that the first coding sequence of a gene starts with a start codon and that the last coding sequence of a gene ends with a stop codon.

## Report your findings in text, and provide the code for the commands you ran in your markdown report.


Command:
```

```
Output:
```

```

Command:
```

```
 Output:
```
 
```
