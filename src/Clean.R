#!/usr/bin/env Rscript

#Function to automate quality control per group with FASTQC/MULTIQC
qualityControl <- function(phenodata,type) {

  # For each group
  for (g in levels(factor(phenodata$Group))) {

    #Output directory for FASTQ reports, one per group
    dir = paste(g,"_qc_",type,sep="")
    cmd = paste("mkdir",dir)
    cat(cmd,"\n")
    system(cmd)

    #For every file in group g
    for (f in phenodata$File[which(phenodata$Group %in% g)]) {

      prefix = ""

      #For every read
      if (type == "clean")
        prefix = "_clean"

      #Single-end reads
      file <- paste(f,prefix,".fastq",sep="")
      if (file.exists(file)) {
        cmd = paste("fastqc -o", dir, file)
        cat(cmd,"\n")
        system(cmd)
      }
      #Else, paired-end reads
      else {
        for (i in 1:2) {
          file <- paste(f,prefix,"_",i,".fastq",sep="")
          cmd = paste("fastqc -o", dir, file)
          cat(cmd,"\n")
          system(cmd)
        }
      }
    }

    #Grouping files with MultiQC
    mqcdir = paste("MQC_",dir,sep="")

    cmd = paste("multiqc -n ",mqcdir," ",dir,sep="")
    cat(cmd,"\n")
    system(cmd)

    cmd = paste("mv ",mqcdir,"_data"," ",dir,sep="")
    cat(cmd,"\n")
    system(cmd)

    cmd = paste("mv ",mqcdir,".html"," " ,dir,sep="")
    cat(cmd,"\n")
    system(cmd)
  }
}

if (!require("readr")) install.packages("readr")
library(readr)

# phenodata.txt file must contain two columns separated by tabs
# - the first one with the name of samples
# - the second with the group to which each samplebelongs

phenodata <- read_delim("phenodata.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#Quality control of raw data
qualityControl(phenodata,"raw")

# Clean with fastp
for(f in phenodata$File) {

  file <- paste(f,".fastq",sep="")

  #For single-end reads
  if (file.exists(file)) {
    cmd = paste("fastp -c -i ",f,".fastq",
                " -o ", f, "_clean.fastq",
                " -j ", f, "_clean.fastp.json ",
                " -h ", f, "_clean.fastp.html ",
                sep="")
  }
  #Else: paired-end reads
  else {
    cmd = paste("fastp -c -i ",f,"_1.fastq"," -I ",f,"_2.fastq",
              " -o ", f, "_clean_1.fastq"," -O ",f,"_clean_2.fastq",
              " -j ", f, "_clean.fastp.json ",
              " -h ", f, "_clean.fastp.html ",
              sep="")
  }

  cat(cmd,"\n")
  system(cmd)
}

#Move fastp reports to fastp folder
cmd = "mkdir fastp"
cat(cmd,"\n")
system(cmd)

cmd = "mv *.json fastp"
cat(cmd,"\n")
system(cmd)

cmd = "mv *.html fastp"
cat(cmd,"\n")
system(cmd)

#Quality control of clean data
qualityControl(phenodata,"clean")
