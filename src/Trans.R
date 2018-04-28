#!/usr/bin/env Rscript

if (!require("readr")) install.packages("readr")
library(readr)

# phenodata.txt file must contain two columns separated by tabs
# - the first one with the name of sra files that we are going to transform to FASTQ (without extension)
# - the second with the group to which each file  belongs

phenodata <- read_delim("phenodata.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
files <- phenodata$File

for (i in 1:length(files))
  files[i] <- paste(files[i],".sra",sep="")

# Conversion of SRA files to FASTQ using the fastq-dump tool

for(f in files) {
  cmd = paste("fastq-dump --split-3", f)
  cat(cmd,"\n")
  system(cmd)
}
