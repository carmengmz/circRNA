#!/usr/bin/env Rscript

if (!require("readr")) install.packages("readr")
library(readr)

# phenodata.txt file must contain two columns separated by tabs
# - the first one with the SRAN RUN identifier that we are going to download from GEO datasets
# - the second with the group to which each RUN belongs

phenodata <- read_delim("phenodata.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# Build the urls for download and store it in a list
files <- phenodata$File
prefix <- "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/"
lurl <- c()

for (i in files) {
  url <- paste(prefix,substr(i,1,6),"/",i,"/",i,".sra",sep="")
  lurl <- c(lurl,url)
}

# Download the SRAs from the list of urls generated previously
if (!require("HelpersMG")) install.packages("HelpersMG")
library(HelpersMG)

wget(lurl) 
