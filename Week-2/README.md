# Assigment 2

## 1: Tell us a bit about the organism.
The organism analyzed was Rhinopithecus bieti, also known as the black-and-white snub-nosed monkey. It is an endangered primate that primarly lives in high-altitude forests of China's and Tibet. As the name suggests, it is characterized by its snub nose with concave nostril. (I thought its nose was silly and that was why I chose it! :))

Command:
```
$ curl -O ftp://ftp.ensembl.org/pub/current_gff3/rhinopithecus_bieti/Rhinopithecus_bieti.ASM169854v1.115.gff3.gz
$ gunzip Rhinopithecus_bieti.ASM169854v1.115.gff3.gz
```
Output:
```
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
  0     0    0     0    0     0      0      0 --:--:-- --:--:--  0     0    0     0    0     0      0      0 --:--:-- --:--:--  2 15.0M    2  451k    0     0   247k      0  0:01:02  0:00:01 13 15.0M   13 2053k    0     0   723k      0  0:00:21  0:00:02 28 15.0M   28 4345k    0     0  1125k      0  0:00:13  0:00:03 37 15.0M   37 5792k    0     0  1205k      0  0:00:12  0:00:04 46 15.0M   46 7190k    0     0  1225k      0  0:00:12  0:00:05 56 15.0M   56 8604k    0     0  1264k      0  0:00:12  0:00:06 67 15.0M   67 10.1M    0     0  1333k      0  0:00:11  0:00:07 79 15.0M   79 11.9M    0     0  1383k      0  0:00:11  0:00:08 87 15.0M   87 13.1M    0     0  1368k      0  0:00:11  0:00:09 98 15.0M   98 14.8M    0     0  1403k      0  0:00:10  0:00:10100 15.0M  100 15.0M    0     0  1400k      0  0:00:10  0:00:10 --:--:-- 1623k
```

## 2: How many sequence regions (chromosomes) does the file contain? Does that match with the expectation for this organism?
This will be discussed in more detait but the assemble for this organism is not complete and therefore is primary composed of contigs and scaffolds. There are 105,032 sequence reqions, but previous studied suggest that the black-and-white snub-nosed monkey has 44 chromosomes. 

Command:
```
$ grep -v '^#' Rhinopithecus_bieti.ASM169854v1.115.gff3 | cut -f1 | sort | uniq | wc -l
```
Output:
```
  105032
```

## 3: How many features does the file contain?
This file contains 1,208,297 features. This a large amount of features but is not unresonable for a mammalian genome annotation. 

Command:
```
$ grep -v '^#' Rhinopithecus_bieti.ASM169854v1.115.gff3 | wc -l
```
 Output:
```
 1208297
```

## 4:How many genes are listed for this organism?
There are 20,966 genes listed for this organism. 

Command:
```
$ grep -v '^#' Rhinopithecus_bieti.ASM169854v1.115.gff3 | awk '$3 == "gene"' | wc -l

```
Output:
```
   20966
```

## 5: Is there a feature type that you may have not heard about before? What is the feature and how is it defined? (If there is no such feature, pick a common feature.)
When looking at all the feature types, I found that this file contained something called lnc_RNA. This feature is for long non-coding RNA, which are RNA molecultes that are above 200 nucleotides but do not code for proteins. These long non-coding RNAs play an import role in regulation of translation, metabolism and signalling. 

To find these unique features, the following command was used:
```
$ grep -v '^#' Rhinopithecus_bieti.ASM169854v1.115.gff3 | cut -f3 | sort | uniq -c | sort -nr | head -20
```
Output:
```
470342 exon
442608 CDS
105032 scaffold
63140 biological_region
43588 mRNA
25624 five_prime_UTR
20966 gene
17694 three_prime_UTR
8376 ncRNA_gene
3091 lnc_RNA
1916 transcript
1762 snRNA
1314 miRNA
 867 snoRNA
 640 rRNA
 563 pseudogenic_transcript
 563 pseudogene
 120 V_gene_segment
  47 scRNA
  22 tRNA
```

## 6: What are the top-ten most annotated feature types (column 3) across the genome?
See table below for rank and count of the most annotated features.

| Rank | Feature Type       | Count    |
|------|--------------------|----------|
| 1    | exon               | 470,342  |
| 2    | CDS                | 442,608  |
| 3    | scaffold           | 105,032  |
| 4    | biological_region  | 63,140   |
| 5    | mRNA               | 43,588   |
| 6    | five_prime_UTR     | 25,624   |
| 7    | gene               | 20,966   |
| 8    | three_prime_UTR    | 17,694   |
| 9    | ncRNA_gene         | 8,376    |
| 10   | lnc_RNA            | 3,091    |

Command: 
```
$ grep -v '^#' Rhinopithecus_bieti.ASM169854v1.115.gff3 | cut -f3 | sort | uniq -c | sort -nr | head -10
```
Output:
```
470342 exon
442608 CDS
105032 scaffold
63140 biological_region
43588 mRNA
25624 five_prime_UTR
20966 gene
17694 three_prime_UTR
8376 ncRNA_gene
3091 lnc_RNA
```

## 7: Having analyzed this GFF file, does it seem like a complete and well-annotated organism?
Putting all of this information together, I can say with confidence that while this genome is well-annoted it is not commplete. The large number of scaffolds (approximately 105k) indicates the genome assembly is fragmented and has yet to be fully assembled into chromosomes. In fact, results show that all of the sequencing regions are from scaffolds. However, the high gene count (approximately 21k) and detailed features like UTRs, ncRNAs, and pseudogenes indicate the annotation itself is comprehensive.

## 8: Share any other insights you might note.

This GFF3 file can still be highly valuable for gene expression, functional genomics, and evolutionary analysis, even though it's not based on a chromosome-level assembly. The presence of ncRNA_gene, lnc_RNA, miRNA, snoRNA, snRNA, etc., shows that non-coding RNAs are included in the annotation and there are 43,588 mRNA entries and 20,966 gene entries, meaning that most genes have multiple transcripts (which is typical). 
