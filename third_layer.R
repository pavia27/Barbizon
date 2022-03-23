#!/usr/bin/env Rscript
rm(list = ls(all.names = TRUE))

#load in libraries
library.path <- .libPaths()
suppressWarnings(library("argparse", lib.loc = library.path, warn.conflicts = F))

#load in classified_contigs
a<-commandArgs(trailingOnly = TRUE)
b<-a[1]
columnHeaders<-c("Contig","classification","reason","lineage","lineage_scores","root","Taxa_Rank_1","Taxa_Rank_2","Taxa_Rank_3","Taxa_Rank_4","Taxa_Rank_5","Taxa_Rank_6","Taxa_Rank_7","Taxa_Rank_8","Taxa_Rank_9","Taxa_Rank_10","Taxa_Rank_11")
c<-as.matrix(read.delim(b,skip =1,header = F, col.names=columnHeaders))

#edit columns
c<-c[,-c(2:6)]
is.na(c) <- c==''
c[is.na(c)] <-"Unassigned"
d<-gsub("\\(.*", "",c)

#writing output
write.table(c,b,row.names=FALSE,quote = FALSE,sep = "\t")

