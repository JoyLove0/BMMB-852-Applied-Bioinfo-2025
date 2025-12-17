# Week 14: Perform an RNA-Seq differential gene expression study
---
### Deliverables:
- A README.md file that summarizes your analysis (written like a concise paper)
- A Makefile that automates your workflow
- A design.csv file that specifies your experimental groups

### Tasks:

- Generate PCA and heatmap visualizations of your data
- Identify a set of differentially expressed genes or transcripts
- Perform functional enrichment analysis on your differentially expressed genes

Your project should demonstrate a complete RNA-Seq analysis pipeline, from FASTQ files through functional enrichment. Keep your README concise and focused—aim for clarity over exhaustive detail.

---
# Setup 
## Stats Environment
```sh
conda activate bioinfo
bio code
bash src/setup/init-stats.sh
```

### Double Check Installation
```sh
micromamba run -n stats Rscript src/setup/doctor.r
```

Output:
```sh
# Doctor, Doctor! Give me the R news!
# Checking DESeq2      ... OK
# Checking gplots      ... OK
# Checking biomaRt     ... OK
# Checking tibble      ... OK
# Checking dplyr       ... OK
# Checking tools       ... OK
# You are doing well, Majesty!
 0
```
# Making the design file
```
 make design
```

# Getting the reference genome and reads

I am using the Airway RNA-Seq study for this assigment. For reference, more information about the data is can be found at:
- [RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC4057123/)
- [Link to Data (via Bioconductor)](https://bioconductor.org/packages/release/data/experiment/html/airway.html)

For my analysis, I "frankensteined" my work from Week-13, the airway.mk, and salmon.mk 


```
make get_data get_genome
```

### Output Summary Tables:
```
processed files:  16 / 16 [======================================] ETA: 0s. done
file                      format  type   num_seqs     sum_len  min_len  avg_len  max_len
reads/SRR1039508_1.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039508_2.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039509_1.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039509_2.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039512_1.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039512_2.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039513_1.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039513_2.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039516_1.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039516_2.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039517_1.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039517_2.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039520_1.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039520_2.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039521_1.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
reads/SRR1039521_2.fastq  FASTQ   DNA   1,000,000  63,000,000       63       63       63
```
```
file                                              format  type  num_seqs      sum_len  min_len  avg_len  max_len
refs/hsapiens/Homo_sapiens.GRCh38.cdna.all.fa.gz  FASTA   DNA    205,541  386,406,247        8  1,879.9  109,224
```

# Indexing and Aligning the Reads
```
make index align 
```

# Create the Counts Matrix and PCA Plots
```
make counts
```

Here is an excerpt of the airway-counts.csv. It's too large to keep the whole thing.

| name           | gene   | N61311_Ctrl | N052611_Ctrl | N080611_Ctrl | N061011_Ctrl | N61311_Dex | N052611_Dex | N080611_Dex | N061011_Dex |
|----------------|--------|-------------|--------------|--------------|--------------|------------|-------------|-------------|-------------|
| ENSG00000000003 | TSPAN6 | 22          | 40           | 39           | 24           | 14         | 23          | 47          | 31          |
| ENSG00000000005 | TNMD   | 0           | 0            | 0            | 0            | 0          | 0           | 0           | 0           |
| ENSG00000000419 | DPM1   | 17          | 24           | 22           | 27           | 32         | 30          | 20          | 20          |
| ENSG00000000457 | SCYL3  | 12          | 12           | 8            | 12           | 13         | 10          | 14          | 19          |
| ENSG00000000460 | FIRRM  | 5           | 0            | 7            | 3            | 7          | 3           | 5           | 4           |
| ENSG00000000938 | FGR    | 0           | 0            | 0            | 0            | 0          | 0           | 0           | 0           |
| ENSG00000000971 | CFH    | 136         | 208          | 268          | 222          | 191        | 284         | 343         | 355         |
| ENSG00000001036 | FUCA2  | 76          | 105          | 88           | 97           | 90         | 92          | 68          | 66          |
| ENSG00000001084 | GCLC   | 28          | 22           | 31           | 28           | 17         | 35          | 25          | 28          |
| ENSG00000001167 | NFYA   | 19          | 26           | 26           | 22           | 19         | 10          | 22          | 11          |
| ENSG00000001460 | STPG1  | 8           | 9            | 7            | 15           | 5          | 2           | 10          | 12          |
| ENSG00000001461 | NIPAL3 | 104         | 163          | 108          | 126          | 81         | 164         | 85          | 149         |
| ENSG00000001497 | LAS1L  | 24          | 19           | 28           | 23           | 16         | 21          | 35          | 27          |
| ENSG00000001561 | ENPP4  | 3           | 5            | 0            | 7            | 2          | 6           | 4           | 6           |
| ENSG00000001617 | SEMA3F | 28          | 26           | 23           | 41           | 17         | 32          | 14          | 27          |
| ENSG00000001626 | CFTR   | 0           | 1            | 0            | 0            | 0          | 0           | 0           | 0           |
| ENSG00000001629 | ANKIB1 | 90          | 83           | 91           | 80           | 63         | 65          | 85          | 58          |
| ENSG00000001630 | CYP51A1| 44          | 40           | 29           | 37           | 41         | 39          | 43          | 34          |
| ENSG00000001631 | KRIT1  | 27          | 41           | 31           | 32           | 38         | 24          | 46          | 30          |
| ENSG00000002016 | RAD52  | 10          | 13           | 10           | 8            | 7          | 4           | 6           | 4           |
| ENSG00000002330 | BAD    | 34          | 28           | 30           | 37           | 38         | 28          | 31          | 24          |
| ENSG00000002549 | LAP3   | 84          | 34           | 38           | 32           | 63         | 72          | 54          | 42          |
| ENSG00000002586 | CD99   | 325         | 325          | 421          | 324          | 352        | 395         | 368         | 314         |

See the PCA plots below:

# Creating Differential Expression Table and Heatmap

```
make edger
```
The top 10 genes with the lowest FDR is as follows:
| name           | gene    | FDR   |
|----------------|---------|-------|
| ENSG00000152583 | SPARCL1 | 0     |
| ENSG00000096060 | FKBP5   | 0     |
| ENSG00000178695 | KCTD12  | 1e-04 |
| ENSG00000166741 | NNMT    | 1e-04 |
| ENSG00000189221 | MAOA    | 1e-04 |
| ENSG00000124766 | SOX4    | 2e-04 |
| ENSG00000148175 | STOM    | 2e-04 |
| ENSG00000162614 | NEXN    | 3e-04 |
| ENSG00000101347 | SAMHD1  | 3e-04 |
| ENSG00000172986 | GXYLT2  | 3e-04 |


If I were to set an arbitrarily limit for FDR to 0.05, that means that all of those genes are diffentially expresses. 

See the Heatmap below:


# Functional Enrichment Analysis
```
make enrich
```

These were several genes that gprofiler found but I only keep genes that had to with cell matrixes, vessels, and lungs because that is biological relevant to this paper. See below:

| ID              | Description                                |
|-----------------|--------------------------------------------|
| GO:0030198      | extracellular matrix organization          |
| GO:0007160      | cell-matrix adhesion                       |
| GO:0001952      | regulation of cell-matrix adhesion         |
| GO:0001568      | blood vessel development                   |
| GO:0048514      | blood vessel morphogenesis                 |
| GO:0031012      | extracellular matrix                       |
| GO:0062023      | collagen-containing extracellular matrix   |
| GO:0005201      | extracellular matrix structural constituent|
| GO:0050840      | extracellular matrix binding               |
| HPA:0300000     | Lung                                       |
| HPA:0301361     | Lung; alveolar cells type I[≥Low]          |
| HPA:0300201     | Lung; endothelial cells[≥Low]              |
| HPA:0301371     | Lung; alveolar cells type II[≥Low]         |
| HPA:0301362     | Lung; alveolar cells type I[≥Medium]       |
| KEGG:05222      | Small cell lung cancer                     |
| KEGG:05223      | Non-small cell lung cancer                 |
| REAC:R-HSA-1474244 | Extracellular matrix organization       |
| WP:WP4658       | Small cell lung cancer                     |


