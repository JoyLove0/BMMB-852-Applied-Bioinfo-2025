# Week 13: Generate an RNA-Seq count matrix
# Tasks: 
- You will need a genome reference and a GTF/GFF annotation file.
- Select at least 6 SRR datasets from the same project, 3 for control and 3 for treatment.
- Write the Makefile that aligns the reads to the genome and creates BAM and BigWig files.
- Run a feature counter to create a count matrix for your data. The final result of your code should be a count matrix that summarizes read counts for each dataset.
- Include IGV screenshots that demonstrate your data is RNA-Seq data.
- Discuss a few lines of the resulting count matrix. Visually identify rows where the counts show consistent gene expression levels. Try to visually verify that the counts in the count matrix are consistent with the numbers you can observe in the alignment tracks.
- Needs to include a README.md, a Makefile, and a design file. Produce a count matrix and IGV screenshots that support the values observed in the count matrix.

# Setup 
## Stats Environment
```sh
bio code
bash src/setup/init-stats.sh
```
# Making the design file
```
 make design
```

Taking a peak with:
```
cat design.csv | column -t -s ,
```

Output:
```
sample  group
HBR_1   HBR
HBR_2   HBR
HBR_3   HBR
UHR_1   UHR
UHR_2   UHR
UHR_3   UHR
```
# Getting the genome reference and a GTF/GFF annotation file
The data is can be found at:
- [Informatics for RNA-seq: A web resource for analysis on the cloud PLoS Computational Biology (2015) by Malachi Griffith, Jason R. Walker, Nicholas C. Spies, Benjamin J. Ainscough, Obi L. Griffith](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004393)
- [Subset of Data](http://data.biostarhandbook.com/data/uhr-hbr.tar.gz)

```
make get_data
```

Output Summary Tables:
```
processed files:  2 / 2 [======================================] ETA: 0s. done
file                       format  type  num_seqs     sum_len     min_len     avg_len     max_len
refs/chr22.genome.fa       FASTA   DNA          1  50,818,468  50,818,468  50,818,468  50,818,468
refs/chr22.transcripts.fa  FASTA   DNA      4,506   7,079,970          33     1,571.2      84,332
```
```
processed files:  6 / 6 [======================================] ETA: 0s. done
file               format  type  num_seqs     sum_len  min_len  avg_len  max_len
reads/HBR_1_R1.fq  FASTQ   DNA    118,571  11,857,100      100      100      100
reads/HBR_2_R1.fq  FASTQ   DNA    144,826  14,482,600      100      100      100
reads/HBR_3_R1.fq  FASTQ   DNA    129,786  12,978,600      100      100      100
reads/UHR_1_R1.fq  FASTQ   DNA    227,392  22,739,200      100      100      100
reads/UHR_2_R1.fq  FASTQ   DNA    162,373  16,237,300      100      100      100
reads/UHR_3_R1.fq  FASTQ   DNA    185,442  18,544,200      100      100      100
```

# Aligning the reads
```
make index align 
```

# Visualizing the alignments with IgV
Bigwig files of all samples:
<img width="1424" height="427" alt="Screenshot 2025-12-11 at 12 55 18 AM" src="https://github.com/user-attachments/assets/ed5bd954-7c84-407f-aae3-796b07babd5b" />

Just HBR_1 and zoomed in:
<img width="1440" height="597" alt="Screenshot 2025-12-11 at 12 57 52 AM" src="https://github.com/user-attachments/assets/86a0f4fc-26a7-445e-8a3a-a078dafc3565" />


# Create the counts
```
make counts
```

Count Summary:
| Status                          | bam/HBR_1.bam | bam/HBR_2.bam | bam/HBR_3.bam | bam/UHR_1.bam | bam/UHR_2.bam | bam/UHR_3.bam |
|---------------------------------|---------------|---------------|---------------|---------------|---------------|---------------|
| Assigned                        | 38511         | 47476         | 41960         | 66216         | 43857         | 52836         |
| Unassigned_Unmapped             | 56043         | 68796         | 61421         | 118463        | 79157         | 98418         |
| Unassigned_Read_Type            | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_Singleton            | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_MappingQuality       | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_Chimera              | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_FragmentLength       | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_Duplicate            | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_MultiMapping         | 1629          | 1968          | 1800          | 7026          | 5478          | 5631          |
| Unassigned_Secondary            | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_NonSplit             | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_NoFeatures           | 19759         | 23258         | 21614         | 33279         | 32899         | 26501         |
| Unassigned_Overlapping_Length   | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_Ambiguity            | 3598          | 4505          | 4069          | 6147          | 3994          | 5056          |


Strand Specific Summary:
| Status                         | bam/HBR_1.bam | bam/HBR_2.bam | bam/HBR_3.bam | bam/UHR_1.bam | bam/UHR_2.bam | bam/UHR_3.bam |
|--------------------------------|---------------|---------------|---------------|---------------|---------------|---------------|
| Assigned                       | 39287         | 48402         | 42852         | 67211         | 43801         | 53605         |
| Unassigned_Unmapped            | 56043         | 68796         | 61421         | 118463        | 79157         | 98418         |
| Unassigned_Read_Type           | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_Singleton           | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_MappingQuality      | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_Chimera             | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_FragmentLength      | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_Duplicate           | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_MultiMapping        | 1629          | 1968          | 1800          | 7026          | 5478          | 5631          |
| Unassigned_Secondary           | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_NonSplit            | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_NoFeatures          | 20796         | 24604         | 22764         | 34566         | 34459         | 27544         |
| Unassigned_Overlapping_Length  | 0             | 0             | 0             | 0             | 0             | 0             |
| Unassigned_Ambiguity           | 1785          | 2233          | 2027          | 3865          | 2490          | 3244          |


I didn't include the formated table because it is too big. Here is a snippet

| name               | gene              | HBR_1 | HBR_2 | HBR_3 | UHR_1 | UHR_2 | UHR_3 |
|--------------------|-------------------|-------|-------|-------|-------|-------|-------|
| ENSG00000100281.13 | ENSG00000100281   | 31    | 45    | 48    | 122   | 99    | 101   |
| ENSG00000273176.1  | ENSG00000273176   | 2     | 0     | 1     | 0     | 0     | 0     |
| ENSG00000100284.20 | ENSG00000100284   | 98    | 119   | 104   | 84    | 55    | 68    |
| ENSG00000266320.1  | ENSG00000266320   | 1     | 0     | 0     | 1     | 1     | 0     |
| ENSG00000274280.1  | ENSG00000274280   | 0     | 0     | 0     | 0     | 0     | 0     |
| ENSG00000282602.1  | ENSG00000282602   | 0     | 0     | 0     | 1     | 0     | 0     |
| ENSG00000282041.1  | ENSG00000282041   | 0     | 0     | 0     | 0     | 0     | 0     |
| ENSG00000100292.16 | ENSG00000100292   | 13    | 18    | 13    | 82    | 51    | 81    |
| ENSG00000100297.15 | ENSG00000100297   | 28    | 31    | 34    | 537   | 292   | 443   |
| ENSG00000233388.2  | ENSG00000233388   | 0     | 0     | 0     | 4     | 1     | 5     |
| ENSG00000273044.1  | ENSG00000273044   | 0     | 0     | 0     | 0     | 0     | 0     |
| ENSG00000100302.6  | ENSG00000100302   | 55    | 81    | 60    | 15    | 14    | 17    |
| ENSG00000198125.12 | ENSG00000198125   | 0     | 1     | 0     | 9     | 8     | 8     |
| ENSG00000268818.2  | ENSG00000268818   | 0     | 0     | 0     | 0     | 0     | 0     |
| ENSG00000221963.5  | ENSG00000221963   | 64    | 54    | 56    | 331   | 197   | 254   |
| ENSG00000228587.1  | ENSG00000228587   | 0     | 0     | 0     | 0     | 0     | 0     |
| ENSG00000226336.2  | ENSG00000226336   | 0     | 0     | 0     | 0     | 0     | 0     |
| ENSG00000128313.2  | ENSG00000128313   | 0     | 0     | 1     | 0     | 0     | 0     |
| ENSG00000100320.22 | ENSG00000100320   | 615   | 733   | 699   | 517   | 362   | 407   |


Several genes show consistent expression within each condition, indicating good replicate agreement. For example, ENSG00000100284 and ENSG00000100320 have similar counts across HBR replicates and across UHR replicates, suggesting stable expression. Some genes display clear differences between HBR and UHR, such as ENSG00000100297 and ENSG00000221963, which are lowly expressed in HBR but strongly expressed in UHR, consistent with condition-specific upregulation. Conversely, ENSG00000100302 shows higher expression in HBR than in UHR, suggesting downregulation in UHR.

Many genes have very low or zero counts across all samples (e.g., ENSG00000274280, ENSG00000273044), indicating little to no expression in either condition; these would typically be filtered out before differential expression analysis. Overall, the table shows a mixture of stably expressed genes, condition-dependent genes, and unexpressed genes, which is expected in an RNA-seq experiment.

I used ENSG00000221963 to vizualize the different expression in IgV. See below:

<img width="1435" height="689" alt="Screenshot 2025-12-16 at 8 22 44 PM" src="https://github.com/user-attachments/assets/9647bf5d-0085-4b23-b7c0-cdc64c58fb87" />



