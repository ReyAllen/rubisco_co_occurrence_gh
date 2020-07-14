# the purpose of this script is to take the summary_results.txt file that is oxyphen-mulitome output,
# read it into dataframe, clean it, and then merge with Form OrgID files in order to correlate
# phylogenetic info with aerobe / anaerobe status.

library(phylotools)
library(tidyverse)
library(dbplyr)
library(tidyr)

setwd("/Users/tonglen/Downloads/oxyphen-multiome/OUTPUT")
oxygen <- as.data.frame(read.table("results_table.tsv", sep = '\t'))

names(oxygen) <- c("filenom", "oxyphen_tot", "oxyphen_ECs")
names(oxygen)
head(oxygen)
oxygen <- separate(oxygen, filenom, c("org_id", "fasta"), sep = "[.]")
head(oxygen)
oxygen <- oxygen[,c(1,3:4)]
dim(oxygen)
head(oxygen)


setwd("/Users/tonglen/Downloads/oxyphen-multiome/Org_ID_inputs/")
IIIb <- read.csv("prok_counts_screened_unique_FormIIIb.csv", sep = "\t")
names(IIIb)
head(IIIb)
head(oxygen)
oxyphen <- merge(oxygen, IIIb, by.x = "org_id", by.y = "Organism_ID")
oxyphen <- oxyphen[,c(2:15,1)]  #move list of EC to end of dataframe because they take up multiple columns when viewed in Open Office
head(oxyphen)
setwd("/Users/tonglen/Downloads/oxyphen-multiome/oxyphen_phylogeny")
write.csv(oxyphen, "FormIIIa_oxyphen.csv")

#////////////////////  code for finding diagnostic EC hits within results using grepl()

