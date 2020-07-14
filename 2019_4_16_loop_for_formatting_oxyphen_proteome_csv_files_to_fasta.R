#purpose of this script is to iterate through a folder of proteome download .csv files, 
# and format them to .fasta files to use as input for OxyPhen


library(phylotools)
library(tidyr)
library(dplyr)

setwd("/Users/tonglen/Downloads/oxyphen-multiome/proteome_downloads/")
filenames <- list.files(path="/Users/tonglen/Downloads/oxyphen-multiome/proteome_downloads/", pattern="*.csv", full.names=FALSE, recursive=FALSE)
head(filenames)
for (i in 1:length(filenames))
{
  handle <- unlist(strsplit(filenames[i], "\\."))  #split filename into org_id and .csv
  handle <- handle[-2] #remove .csv from vector filenames[i]
  proteome <- read.csv(filenames[i], sep = '\t')
  proteome <- subset(proteome, select=c(Organism, Entry.name, Protein.names, Sequence))
  proteome <- unite(proteome, "Nom", c(Organism, Entry.name, Protein.names), sep="_")
  names(proteome) <- c("seq.name", "seq.text")
  dat2fasta(proteome, outfile = (sprintf("/Users/tonglen/Downloads/oxyphen-multiome/PROTEOMES/%s.fasta", handle)))
}

