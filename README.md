# Classification of RNA-Seq samples using circRNA expression
A collection of R-scripts to automate the download of RNA-Seq samples, quality control, preprocessing, detection of circRNAs and classification using machine learning from the detected circRNAs.

## User manual
In first place we must select a set of sequenced RNA-seq samples suitable for the detection of circRNAs. In [Selection of the RNA-Seq library](https://github.com/carmengmz/circRNA/wiki/Selection-of-the-RNA-seq-library) you can find some tips on how to make the selection.

If the detected circRNAs will be used for classification, the samples must belong to different groups. For example, sequenced samples from normal patients vs sequenced samples from patients with some disease or condition. In our example we will select sequenced samples of blood exosomes from normal patiens and sequenced samples of blood exosomes from patiens with coronary heart disease.

We have performed the search in the [Gene Expression Omnibus repository (GEO)](https://www.ncbi.nlm.nih.gov/gds) of the NCBI (National Center for Biotechnology Information) using the following criteria:
- <b>Organism</b>: Homo sapiens
- <b>Study type</b>: Expression profiling by throughput sequencing
- <b>Query</b>: "blood exosomes"

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

Now we can run the script [Download.R](https://github.com/carmengmz/circRNA/blob/master/src/Download.R) to download the <b>*.sra</b> files. The <b>phenodata.txt</b> file must be in the working directory. In a Unix like S.O. command line we will type:
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
Once the SRAs are downloaded, the next step is to convert the files to the FASTQ format, separating the readings into two different files in the case of paired-end reads. To do this we are going to use the <b>fastq-dump</b> tool included in the [NCBI's SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) available for several operating systems. The automation of this process has been implemented in the [Trans.R](https://github.com/carmengmz/circRNA/blob/master/src/Trans.R) script. The command to convert the files to the FASTQ format with <b>fastq-dump</b> is:
```
fastq-dump --split-3 <file.sra>
```
The <b>--split-3</b> option generates a single FASTQ file if we have single-end reads or two FASTQ files separating the readings paired in case of paired-end reads. In this case the file names will have the suffixes _1.FASTQ and _2.FASTQ. In addition, we can obtain a third file FASTQ without suffix and of a much smaller size, which will contain "orphan" readings. We will discard this file.

The <b>fastq-dump</b> tool must be able to be executed directly in the working directory, for example, adding the path to the executable to the <b>.bashrc</b> file on a Unix system or to the global PATH variable on a Windows system. The <b>*.sra</b> files and the <b>phenodata.txt</b> file must be too in the working directory. In a Unix like S.O. command line we will type:

```
> Rscript Trans.R
```
