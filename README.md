# Classification of RNA-Seq samples using circRNA expression
A collection of R-scripts to automate the download of RNA-Seq samples, quality control, preprocessing, detection of circRNAs and classification using machine learning from the detected circRNAs.

## User manual
In first place we must select a set of sequenced RNA-seq samples suitable for the detection of circRNAs. In [Selection of the RNA-Seq library](https://github.com/carmengmz/circRNA/wiki/Selection-of-the-RNA-seq-library) you can find some tips on how to make the selection.

If the detected circRNAs will be used for classification, the samples must belong to different groups. For example, sequenced samples from normal patients vs sequenced samples from patients with some disease or condition. In our example we will select sequenced samples of blood exosomes from normal patiens and sequenced samples of blood exosomes from patiens with coronary heart disease.

We performed the search in the [Gene Expression Omnibus repository (GEO)](https://www.ncbi.nlm.nih.gov/gds) of the NCBI (National Center for Biotechnology Information) using the following criteria:
- Organism: Homo sapiens
- Study type: Expression profiling by throughput sequencing
- Query: ""blood exosomes"

From the search results we have selected two projects containing sequenced samples from human blood exosomes, which were also used to characterize circRNA in its associated publication.
- [RNA-seq reveals abundant circRNA, lncRNA and mRNA in blood exosomes of patients with coronary heart disease](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99985)
- [RNA-seq reveals abundant circRNA, lncRNA and mRNA in blood exosomes of normal persons](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100206)

From every result (follow the links) we retrieve the id of the BioProyect: PRJNA390278 and 	PRJNA390988.

Once we have the BioProyect id, we will retrieve the ids of the SRA Runs in: https://www.ncbi.nlm.nih.gov/Traces/study/ (under the column Run) and we will build a text file separated by tabs called <b>phenodata.txt</b> with two columns: Name and Group. In our [example](https://github.com/carmengmz/circRNA/tree/master/example) we will work with three samples of each group:

```
File  Group
SRR5679909  coronary
SRR5679908  coronary
SRR5679907  coronary
SRR5712482  normal
SRR5712483  normal
SRR5712484  normal
```

Now we can run the script [Download.R](https://github.com/carmengmz/circRNA/blob/master/src/Download.R) to download the *.sra files. The phenodata.txt file must be in the working directory. In a Unix like S.O. command line we will type:
```
> Rscript Download.R
```
As a result, at the end of the script, in the working directory we will have the downloaded sra files and a log file:
```
download.log
phenodata.txt
SRR5679909.sra
SRR5679908.sra
SRR5679907.sra
SRR5712482.sra
SRR5712483.sra
SRR5712484.sra 
```
Now we can move to the next step.
