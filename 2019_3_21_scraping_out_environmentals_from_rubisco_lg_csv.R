# the purpose of this script is to take the uniprot lg_rub download csv and filter out
# uncultured / environmental samples from the 146,000 uniprot entries.
# also scrape out eukaryota to get a final list of rub_lg prokaryote fasta for annotating isoforms via tree2fasta


setwd("/Users/tonglen/Downloads/rub-co-occurrence-master_2/uniprotKB_downloads")
rub_uniprot <- read.csv("rubisco_large.csv")
dim(rub_uniprot) #159,116
names(rub_uniprot)

#looking for synechocystis:
synecho <- subset(rub_uniprot, rub_uniprot$Organism.ID == "1111708")
synecho


setwd("/Users/tonglen/large_tree/")

arch <- read.delim("assembly_summary_refseq_archaea.txt", sep="\t", stringsAsFactors = FALSE) #   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the columns in this file.
dim(arch)  #901
arch_compl <- arch[arch$assembly_level == "Complete Genome",]
dim(arch_compl)  #299


bact <- read.delim("assembly_summary_refseq_bacteria.txt", sep="\t", stringsAsFactors = FALSE)
dim(bact)  #165,327
#check for synecho 6803:
synecho <- subset(bact, bact$taxid == "1111708")
synecho
bact_compl <- bact[bact$assembly_level == "Complete Genome",]
dim(bact_compl)  #14,379
write.csv(bact_compl, "bact_compl.csv", row.names = F)

prok_compl <- rbind(bact_compl, arch_compl)
dim(prok_compl) #14,678   #list of proks with complete genomes, with or without rubisco in their genomes.
head(prok_compl)
names(prok_compl)
#next, remove any orgid (lg_rub uniprot) that does not correspond with a complete genome taxid (ncbi refseq).

#don't use merge, you don't know how it works
#merge removes duplicates 
#rub_proks_compl <-  merge(rub_proks, prok_compl, by.x="Organism.ID", by.y="taxid") 
#dim(rub_proks_compl)  #8080


#instead use subset(), and use %in% to create a vector / list you want to use to filter the df.
# https://stackoverflow.com/questions/38850629/subset-a-column-in-data-frame-based-on-another-data-frame-list

#rub_proks_compl <- subset(rub_proks, rub_proks$Organism.ID %in% prok_compl$taxid)  #hahah just kidding this one will remove duplicate rows as well :/
#dim(rub_proks_compl)  #1598

#this one says it doesn't remove duplicate rows?
#https://stackoverflow.com/questions/46656244/r-how-to-subset-duplicate-rows-of-data-frame
rub_proks_compl <- rub_uniprot[match(prok_compl$taxid, rub_uniprot$Organism.ID), ] 
#mydata[match(ids, mydata$id), ]  where match(x, table) uses vector x to select any rows with x in table$column
#match(x, table)  gives all the rows in rub_proks that have a taxid x with complete genome, hopefully
#it does this by inserting rows with 'NA' in all the columns, not helpful.
dim(rub_proks_compl) #14,678
head(rub_proks_compl)
rub_proks_compl <- rub_proks_compl[!is.na(rub_proks_compl$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
#this step should keep all rubisco paralogs, but remove rows from org.IDs not found in prok_compl

dim(rub_proks_compl) #1655
names(rub_proks_compl)
head(rub_proks_compl)  

#unique_rub_proks <- rub_proks_compl[!duplicated(rub_proks_compl$Organism.ID), ] 
#subsets a dataframe of unduplicated org.ids
#dim(unique_rub_proks) #911

####eliminated this step Dec 2019 as I think it would inappropriately remove any 
#instances of multiple rub_lg entries in one orgID, when those should be included,
#as many orgs have multiple copies of different rub_lg.
#nums_9 completed this de-replication pipeline; nums_10 will be run without this de-replication.

unique_rub_proks <- rub_proks_compl
#  screen out "uncultured", "environmental", and "metagenomic" orgs from rub_uniprot download

cultured <- subset(unique_rub_proks, !grepl("uncultured", unique_rub_proks$Organism), drop = TRUE)
dim(cultured)  #1655
not_environmental <- subset(cultured, !grepl("environmental", cultured$Organism), drop = T)
dim(not_environmental) #911
not_enrichment <- subset(not_environmental, !grepl("enrichment", not_environmental$Organism), drop = T)
not_metagenomic <- subset(not_enrichment, !grepl("metagenomic", not_enrichment$Organism), drop = T)
not_metagenomic <- subset(not_metagenomic, !grepl("metagenome", not_metagenomic$Organism), drop = T)
dim(not_metagenomic) #1655


#///////////////////////////////  next, remove fragment and short sequences (likely from incomplete genome seqs)
defrag <- subset(not_metagenomic, Fragment!='fragment')  #create defrag df that exludes rubsico seqs that only includes full-length rubisco seq and excludes partial seqs. 
dim(defrag)  #1597

head(defrag)
names(defrag)
library(phylotools)
dat2fasta(defrag, "rub_prok_compl.fasta")  #suprisingly, 70 complete genomes only contain fragment of rubisco protein.
clean.fasta.name("rub_prok_compl.fasta", "rub_prok_compl_cleaned.fasta")
#////////////////////////////  next, remove un-needed ncbi columns


names(defrag) <- c("Entry",	"Entry name",	"Gene names",	"Gene names  (ORF )",	"Organism",	"Organism ID",	"Length",	"Mass",	"Sequence",	"EC number",	"Fragment",	"Pathway",	"Annotation",	"Status",	"Intramembrane",	"Transmembrane",	"Cross reference (INTERPRO)",	"Cross reference (PFAM)",	"Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)",	"Taxonomic lineage (PHYLUM)",	"Taxonomic lineage (CLASS)",	"Taxonomic lineage (ORDER)",	"Cofactor",	"Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							
head(defrag)
write.csv(defrag, "screen_for_length.csv", row.names=FALSE)  #checked in open office; seems fine
head(defrag)

# look for missing synechocystis 6803 orgId in defrag df:

synecho <- subset(defrag, defrag$`Organism ID` == "1111708")
synecho


# ok the two columns we want for annotation are "organism" and "sequence"
library(phylotools)
library(dplyr)
# subset(dat, select=c("A", "B"))
names(defrag)
defrag <- defrag[, c(5,9)] #confirm these are the correct columns to excerpt for fasta
dim(defrag)
#/////////////////////////////////  change col names to convert df to fasta file
names(defrag) <- c("seq.name", "seq.text")
head(defrag)

#merge defrag large_rub with annotated forms for a master defragged list of rubisco-having genomes

setwd("/Users/tonglen/tree/data")
getwd()
library(phylotools)
merged_all <- read.fasta("merged_all.fasta")  #the previous fasta file used for pipeline, includes euks and uncultureds filtered out later
dim(merged_all)  #920 rows, 2 cols


setwd("/Users/tonglen/large_tree/Tree2Fasta_Output_nums_9/") 

library(phylotools)


annotated_forms <- read.fasta("shih_badger_jaffe_kono_annotated_forms.fas")
head(annotated_forms)

dim(defrag)  #1596
dim(annotated_forms)  #382


nums <- as.data.frame(rnorm(1597, 100, 5))
colnames(nums) <- "Nums"
dim(nums)
head(nums)
cultured_proks_nums <- cbind(nums, defrag)
head(cultured_proks_nums)
library(dplyr)
library(tidyverse)
library(tidyr)
cultured_proks_nums <- unite(cultured_proks_nums, seq.name, c(seq.name, Nums), sep="_")

head(cultured_proks_nums)
merged_all_10 <- rbind(annotated_forms, cultured_proks_nums)
dim(merged_all_10) #1979
head(merged_all_10)

#//////////////////////////////  make fasta file and then clean fasta headers

setwd("/Users/tonglen/large_tree/Tree2Fasta_Output_nums_10/")
dat2fasta(merged_all_10, outfile = "merged_all_10.fasta")
clean.fasta.name(infile = "merged_all_10.fasta", outfile = "merged_all_10_cleaned.fasta")
dim(merged_all_10) #1979


#////////////////////  add random numbers to fasta headers so that each header is unique to avoid errors in tree-building

nums <- as.data.frame(rnorm(894, 100, 5))
colnames(nums) <- "Nums"
dim(nums)
head(nums)
merged_cleaned_nums <- cbind(nums, merged_cleaned)
head(merged_cleaned_nums)
library(dplyr)
library(tidyverse)
library(tidyr)
merged_nums_fasta <- unite(merged_cleaned_nums, seq.name, c(seq.name, Nums), sep="_")

head(merged_nums_fasta)
dat2fasta(merged_nums_fasta, outfile = "merged_cleaned_nums.fasta")
clean.fasta.name(infile = "merged_cleaned_nums.fasta", outfile = "merged_cleaned_nums.fasta")
head(merged_nums_fasta)
head(cleaned_jaffe_forms)
merged_all <- rbind(merged_nums_fasta, cleaned_jaffe_forms)
head(merged_all)

dim(merged_nums_fasta)
dim(cleaned_jaffe_forms)
dim(merged_all)

dat2fasta(merged_all, outfile = "merged_all.fasta")
clean.fasta.name(infile = "merged_all.fasta", outfile = "merged_all.fasta")

#/////////////////////////////  troubleshoot grepl() command, which is including some that should be excluded and vice versa, according to 
# open office control+F searches of resulting files

a = c("unculturedspecies", "ecoli", "ecoli", "ecoli")
b = c("uncultured species", "ecoli", "ecoli", "ecoli")
c = c("uncultured", "ecoli", "ecoli", "ecoli")
d = c("metagenomicsample", "ecoli", "ecoli", "ecoli")
e = c("metagenomic sample", "ecoli", "e coli", "e coli")
f = c("metagenomic", "ecoli", "ecoli", "ecoli")
g = c("environmentalsample", "e coli", "e coli", "e coli")
h = c("environmental sample", "e coli", "e coli", "e coli")
i = c("environmental", "e coli", "e coli", "e coli")
j = c("e coli", "e coli", "e coli", "e coli")

df = as.data.frame(rbind(a,b,c,d,e,f,g,h,i,j))
names <- c("col1", "col2", "col3", "col4")
colnames(df) <- names
df

columns <- c("col1", "col2", "col3", "col4")
metag <- c("environmental", "uncultured", "metagenomic")
#nonenv_df <- NULL
#nonenv_df <- grepl.sub(df, metag, columns, keep.found = F)  #this shit is simply not working at all anymore :/
# and when it does "work", it removes rows it shouldn't and leaves rows that should be removed.

#filtered <- df[!(grepl("uncultured", df$col1) & 
#              (grepl("metagenomic", df$col1) |
 #                grepl("environmental", df$col1))),]

#library(DataCombine)
#subset(data, data[[2]] %in% c("Company Name 09", "Company Name"), drop = TRUE) 
#subset(data, grepl("^Company Name", data[[2]]), drop = TRUE)

dim(rub_proks)


#subset(data, data[2] == "Company Name 09" | "Company Name", drop = T)

#non_env <-  as.data.frame(subset(df, !grepl("uncultured", df$col1), drop = T) %in%
#            subset(df, !grepl("environmental", df$col1), drop = T) %in%
#            subset(df, !grepl("metagenomic", df$col1), drop = T))

#non_env   this only produces a set of logicals, not the subsetted dataframe :Z
#//////////////////////


colnames(rub_proks)

colnames(df)


check <- grepl.sub(rub_proks, metag, columns, keep.found = T)
dim(check)
write.csv(check, "check_uncultured.csv")
dim(proks_nonenv)
head(proks_nonenv)  #organism name still not squished.
write.csv(proks_nonenv, "proks_nonenv_3.csv")  #checked with Find function for uncultured, environmental in OpenOffice
