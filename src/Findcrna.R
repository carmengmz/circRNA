#!/usr/bin/env Rscript

#MODIFY THIS VARS WITH DESIRED VALUES
####################################################
num_threads = 10
ref_fasta = "hg38.fa"
ref_gtf = "hg38.gtf"
####################################################

if (!require("readr")) install.packages("readr")
library(readr)

# phenodata.txt file must contain two columns separated by tabs
# - the first one with the name of files the sra we are going to transform to FASTQ (without extension)
# - the second with the group to which each file  belongs

phenodata <- read_delim("phenodata.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
files <- phenodata$File

cmd = paste("bwa index -a bwtsw ", ref_fasta)
cat(cmd,"\n")
system(cmd)

#Align sequences to the reference genome

for(f in phenodata$File) {

  file <- paste(f,"_clean.fastq",sep="")

  #For single-end reads
  if (file.exists(file)) {
    cmd = paste("bwa mem -T 19 -t ",num_threads ," ", ref_fasta, " ", file, "> aln-", f,".sam", sep="")
  }

  #Else: paired-end reads
  else {
    cmd = paste("bwa mem -T 19 -t ",num_threads ," ", ref_fasta, " ", f, "_clean_1.fastq ",f,"_clean_2.fastq > aln-", f,".sam", sep="")
  }

  cat(cmd,"\n")
  system(cmd)
}

#Detect circRNA with CIRI2

for(f in phenodata$File) {

  file <- paste("aln-",f,".sam",sep="")
  cmd = paste("perl CIRI2.pl -T 10 -I ",file," -O  outfile-",f," -F ", ref_fasta," -A ", ref_gtf, "> ciri-",f,".log", sep="")
  cat(cmd,"\n")
  system(cmd)
}
