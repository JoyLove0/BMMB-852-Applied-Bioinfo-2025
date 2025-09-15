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


## From your GFF file, separate the intervals of type "gene" or "transcript" into a different file. Show the commands you used to do this.
## Visualize the simplified GFF in IGV as a separate track. Compare the visualization of the original GFF with the simplified GFF.
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
