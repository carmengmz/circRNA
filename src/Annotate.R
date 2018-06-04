#!/usr/bin/env Rscript

# MODIFY THIS VAR WITH DESIRED RANGE OF ANNOTATION
####################################################
range = 10
####################################################

if (!require("readr")) install.packages("readr")
library(readr)

# phenodata.txt file must contain two columns separated by tabs
# - the first one with the name of samples
# - the second with the group to which each samplebelongs
phenodata <- read_delim("phenodata.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#Empty annotation dataframe
circRNA_ann <- data.frame(circ=character(), chr=character(), circRNA_start = integer(), circRNA_end = integer(), 
                            strand = character(), stringsAsFactors=FALSE)
Files <- phenodata$File

for (f in Files) {

  name <- paste("outfile-",f,sep="")
  circ <- read_delim(name,"\t", escape_double = FALSE, trim_ws = TRUE)

  for (i in 1:nrow(circ)) {

    #Data for row i in sequence f
    start <- circ$circRNA_start[i]
    end <- circ$circRNA_end[i]
    chrom <- circ$chr[i]
    strand <- circ$strand[i]
    reads <- circ$`#junction_reads`[i]

    #In the first row of first sample: initialize data base
    if (nrow(circRNA_ann)==0) {

      circRNA_ann <- data.frame(as.character("circ1"),as.character(chrom), as.integer(start), as.integer(end),
                                as.character(strand), as.integer(0), stringsAsFactors = FALSE)
      names(circRNA_ann) <- c("circ", "chrom", "start", "end", "strand",f)
      rownames(circRNA_ann) <- circRNA_ann$circ
    }

    #In the first row from the second sample to the end: initialize column for this sample
    else if (i==1) {
      circRNA_ann[, f] <- 0
    }

    # Searching a circRNA in the same chrom and strand with start +/- range and end +/- range
    match <- circRNA_ann[  circRNA_ann$start - range <= start & circRNA_ann$start + range >= start &
                           circRNA_ann$end - range  <= end & circRNA_ann$end + range >= end &
                           chrom == circRNA_ann$chrom & strand == circRNA_ann$strand, ]

    #If empty match: add the circRNA to data base
    if (nrow(match) == 0) {
      circId <- paste("circ",nrow(circRNA_ann)+1,sep="")

      circRNA_ann <- rbind.data.frame(circRNA_ann,
                     c(as.character(circId), as.character(chrom), as.integer(start), as.integer(end),
                       as.character(strand),rep(0,ncol(circRNA_ann)-5)))
      rownames(circRNA_ann) <- circRNA_ann$circ

      for (k in 6:ncol(circRNA_ann))
        circRNA_ann[,k] <- as.integer(circRNA_ann[,k])

      circRNA_ann$start <- as.integer(circRNA_ann$start)
      circRNA_ann$end <- as.integer(circRNA_ann$end)

      circRNA_ann[circId,f] <- as.integer(reads)
    }

    #If there was a match: add the counts to the row in this column sample
    else if (nrow(match==1)) {
      circId <- match$circ[1]
      circRNA_ann[circId,f] <- circRNA_ann[circId,f] + as.integer(reads)
    }

    #More than one match: this should not happen. 
    else if (nrow(match>1)) {
      print(match)
      stop()
    }
  }

  print(paste("Finished",f))
}

saveRDS(circRNA_ann,"circ_anotations.rds")
