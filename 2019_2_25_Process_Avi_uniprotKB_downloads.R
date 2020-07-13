#read Avi's .csv scraped from UniProtKB api
uniprot_rubisco_lg <- read.csv("rubisco_large.csv", sep = "\t", stringsAsFactors = FALSE)
library(dplyr)    
rub_lg_prok <- uniprot_rubisco_lg %>% filter(Taxonomic.lineage..SUPERKINGDOM. %in% c("Bacteria", "Archaea"), Status == "reviewed")
dim(rub_lg_prok)  
# gives 182 entries)
rub_lg_rev <- uniprot_rubisco_lg %>% filter(Status == "reviewed")
dim(rub_lg_rev)
# gives 710 entries.
dim(rubisco_lg)
#146,000 entries
#note that there are some "rubisco-like proteins" in the mix 
write.csv(rub_lg_rev, "rub_lg_rev.csv")  #export filtered file
#rub_lg_rev is the name of the uniprot rubiscos; three columns

#///////////////////////////////
#create fasta files to upload to FigTree etc
# want format: >formI Synechocystis MASLDKFJSLDFKJSDLKJSD
# want format: >EntryNumber Synechocystis MSDLKFJSDLKFJ

all_forms <- read.csv("all_forms_formatted.csv", sep = "\t")
dim(all_forms)
head(all_forms)
library(dplyr)
fasta_forms <- select(all_forms, "FORM", "Organism", "Sequence")
head(fasta_forms)
dim(fasta_forms)
#ok got 275 various formses.
library(tidyr)
headers <- unite(fasta_forms, header, "FORM", "Organism", sep = "_", remove = TRUE) # merge FORM and Organism into single column, remove separate columns
dim(headers)
head(headers)

library(phylotools)

names(headers) <- c("seq.name", "seq.text") #The column of the data frame must be: 1. seq.name, 2. seq.text, represent the name of the sequences, the content of the sequence, eg. ATCGGGAAC. 
#headers is the name of the dataframe with rubisco forms

# //////////////////////////////////
#format rub_lg_rev, UniprotKB excerpt of ~700 reviewed entries

rub_lg <- select(rub_lg_rev, "Organism", "Entry", "Sequence")
head(rub_lg)
dim(rub_lg)
#ok got 710 Uniprot rubiscos.
library(tidyr)
uniprot <- unite(rub_lg, header, "Entry", "Organism", sep = "_", remove = TRUE) # merge Entry and Organism into single column, remove separate columns
head(uniprot)
dim(uniprot)
library(phylotools)

names(uniprot) <- c("seq.name", "seq.text") #The column of the data frame must be: 1. seq.name, 2. seq.text, represent the name of the sequences, the content of the sequence, eg. ATCGGGAAC. 


dat2fasta(uniprot, outfile = "rub_lg.fasta")  #make it into fasta file. 
clean.fasta.name(infile = "out.fasta", outfile = "rubisco_uniprot_cleaned.fasta") #clean fasta header to remove whitespace etc

#////////////////////////////////////////
#mergity-merge the forms dataframe and the uniprot dataframe, and export as fasta


merged_forms_uniprot <- rbind(headers, uniprot, by = c("seq.name", "seq.text"))
dim(merged_forms_uniprot)

dat2fasta(merged_forms_uniprot, outfile = "merged_forms_uniprot.fasta")  #make it into fasta file. 
clean.fasta.name(infile = "merged_forms_uniprot.fasta", outfile = "rubisco_merged_cleaned.fasta") #clean fasta header to remove whitespace etc

#///////////////////////////////////////
