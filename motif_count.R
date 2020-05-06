#!/usr/bin/env Rscript
# Use: Rscript motif_count.R input.fa motif max.mismatches output.txt

suppressMessages(library(seqinr))
suppressMessages(library(Biostrings))

args = commandArgs(trailingOnly=TRUE)

dna <- read.fasta(args[1], seqonly = TRUE)
motif <- args[2]
match <- args[3]
output <- args[4]

df <- NULL
for (i in 1:length(dna)) { 
      x <- dna[[i]];
      count <- countPattern(motif, x, max.mismatch=as.integer(match), min.mismatch=0, with.indels=TRUE, fixed=TRUE, algorithm="auto");
      df <- rbind(df, data.frame(x,count))
    }

write.table(df, file=output, quote=FALSE, sep="\t", row.names=FALSE)

