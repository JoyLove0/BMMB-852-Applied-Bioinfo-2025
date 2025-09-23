# Assignment 4: Obtain genomic data via accession numbers

I am in group 2 and used the paper titled: "Zika Virus Targets Human Cortical Neural Precursors and Attenuates Their Growth". This paper can be found at this link: https://pubmed.ncbi.nlm.nih.gov/26952870/.

## Identify the accession numbers for the genome and write the commands to download the data. Make your commands reusable and reproducible.
For this assignment, the Zika genome is needed.

Command: 
```
datasets download genome accession GCF_000882815.3 --include genome,gff3,gtf
unzip ncbi_dataset.zip
cd ncbi_dataset/data/GCF_000882815.3/
mv GCF_000882815.3_ViralProj36615_genomic.fna zika_paper_genome.fa
```

## Now use IGV to visualize the genome and the annotations relative to the genome.

Genome:
Genomes > Load Genome From File > Selected zika_paper_genome.fa

Annotation Files:
File > Load from file > Selected genomic.gtf 
File > Load from file > Selected genomic.gff 

<img width="1440" height="900" alt="Screenshot 2025-09-23 at 4 59 57 PM" src="https://github.com/user-attachments/assets/86a6c34e-9e59-4318-8b6e-72ff57fde89c" />

## How big is the genome? 
The genome is 10794 bps. 
Command:
```
awk '!/^>/ {sum += length($0)} END {print sum}' zika_paper_genome.fa
```
Output
```
10794
```

## How many features of each type does the GFF file contain? 
| Count | Feature                    |
|-------|----------------------------|
| 14    | mature_protein_region_of_CDS |
| 1     | three_prime_UTR           |
| 1     | region                    |
| 1     | gene                      |
| 1     | five_prime_UTR           |
| 1     | CDS                       |

Command:
```
awk '$0 !~ /^#/ {print $3}' genomic.gff | sort | uniq -c | sort -nr
```
Output
```
14 mature_protein_region_of_CDS
   1 three_prime_UTR
   1 region
   1 gene
   1 five_prime_UTR
   1 CDS
```

## What is the longest gene? What is its name and function? You may need to search other resources to answer this question. 
The longest gene is gene-ZIKV_gp1 and it encodes for the GP1 protein (an envolpe protein) that enables the virus to bind to host receptors. 

Command:
```
awk '$3 == "gene" && $0 !~ /^#/ && match($9, /ID=([^;]+)/, a) { print a[1], $5 - $4 + 1 }' genomic.gff
```

Output
```
gene-ZIKV_gp1 10260
```

## Look at other gene names. Pick another gene and describe its name and function.
There is only one gene. I checked with the following command:
```
awk '$0 !~ /^#/ && $3 == "gene"' genomic.gff 
```

Output:
```
NC_012532.1	RefSeq	gene	107	10366	.	+	.	ID=gene-ZIKV_gp1;Dbxref=GeneID:7751225;Name=POLY;gbkey=Gene;gene=POLY;gene_biotype=protein_coding;locus_tag=ZIKV_gp1
```
## Look at the genomic features, are these closely packed, is there a lot of intragenomic space? 

This genome is not packed at all- considering there is only 1 gene and all features are associated with this gene.

## Using IGV estimate how much of the genome is covered by coding sequences.

Almost the whole genome is covered by the coding region:
<img width="1440" height="900" alt="Screenshot 2025-09-23 at 5 54 34 PM" src="https://github.com/user-attachments/assets/2dad7f56-4a4f-4501-ba1f-b04ea0a460ce" />

## Find alternative genome builds that could be used to perhaps answer a different question (find their accession numbers). Considering the focus of the paper, think about what other questions you could answer if you used a different genome build.

I went to https://www.ncbi.nlm.nih.gov/datasets/ and searched for alternative Zika genomes. There are only five genomes avaliable there and most them are associated with other strains of Zika. 

<img width="1430" height="802" alt="Screenshot 2025-09-23 at 5 59 55 PM" src="https://github.com/user-attachments/assets/e847c336-4c4f-4440-bcde-53f1e7f5b466" />

Since there are several strains of this same virsuses, comparitive strain studies could be conducted. This would allow us to answer questions like: How do different zika strains compare in their ability to infect human cortical neural progenitor cells? What genetic variations exist between your strains and how do those mutations effect virulence or neurotropism?

