# purpose of this script is to take the sorted rubsico seqs that are Tree2Fasta output by FormIAc, FormIAq, etc.
# use that to generate uniprot files of orgi_id, phylogeny for each Form.
# output of this script will be inputs for avi's python code for calculating co-occurrence and heat map by rubisco form.
# include screening by the ncbi refseq list of taxid with completely sequenced genomes.
# output of this script also used as input for oxyphen analysis.

#  "Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway	Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names"																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							


# read in Tree2Fasta output files.

rm(list=ls())

setwd("/Users/tonglen/large_tree/Tree2Fasta_Output_nums_10/FASTA_by_annot_")
library(phylotools)
library(dplyr)
library(tidyverse)

FormI <- read.fasta("I.fas")
FormIA <- read.fasta("IA.fas")
FormIAc <- read.fasta("IAc.fas")
FormIAq <- read.fasta("IAq.fas")
FormIB <- read.fasta("IB.fas")
FormIBc <- read.fasta("IBc.fas")
FormIC <-  read.fasta("IC.fas")
FormIE <- read.fasta("IE.fas")
FormII <- read.fasta("II.fas")
II_III <- read.fasta("II_III.fas")
FormIII <- read.fasta("III.fas")
FormIIIa <- read.fasta("IIIa.fas")
FormIIIb <- read.fasta("III_b.fas")
FormIIIc <- read.fasta("III_c.fas")
FormIV <- read.fasta("IV.fas")



setwd("/Users/tonglen/large_tree/")
library(phylotools)

proks_compl <- read.csv("screen_for_length.csv")  #defragged complete genomes using refseq prok files.

dim(proks_compl)   #1597
head(proks_compl)
names(proks_compl) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

#names(proks_compl) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway", "Annotation", "Status", "Intramembrane", "Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)", "Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")
head(proks_compl)
setwd("/Users/tonglen/Downloads/rub-co-occurrence-master_2/uniprotKB_downloads/isoforms") #move to folder for input files for co-occurrence analysis

# merge() function: "The rows in the two data frames that match on the specified columns are extracted, 
# and joined together. If there is more than one match, all possible matches contribute one row each."
#it *sounds* like duplicates are not removed? But I can't really tell. 



FormI_uni <- merge(proks_compl, FormI, by.x = "Sequence", by.y="seq.text", all.y=TRUE) #merge FormI fasta with download csv by sequence, does not remove identical sequences (paralogs)
dim(FormI_uni)  #19 26
head(FormI_uni) 
head(FormI)  #fasta file with seq.name and seq.txt as columns.
FormI_uni_uni <- FormI_uni[ ,c(2:9, 1, 10:25)]  #re-arrange columns to put them in same order as uniprot downloads. 
names(FormI_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							
dim(FormI_uni_uni) #9 25
FormI_uni_uni <-  FormI_uni_uni[!duplicated(FormI_uni_uni[,'Entry']),]  #removes duplicate rows as determined by the "Entry" column. Does not remove paralogs.
dim(FormI_uni_uni)  #12 25
FormI_uni_uni <- FormI_uni_uni[!is.na(FormI_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormI_uni_uni) #12 25

head(FormI_uni_uni)  

write.csv(FormI_uni_uni, "FormI.csv", row.names=FALSE) 

dim(FormIA)

FormIA_uni <- merge(proks_compl, FormIA, by.x = "Sequence", by.y="seq.text", all.y=TRUE) 
dim(FormIA_uni)  #12 
head(FormIA_uni)
FormIA_uni_uni <-  FormIA_uni[!duplicated(FormIA_uni[,'Entry']),]  #removes duplicate rows as determined by the "Entry" column. Does not remove paralogs.
dim(FormIA_uni_uni) #4  26  should be smaller than FormI, because there are duplicate seqs in the (pre-annotated and post-annotated) fasta files. All were appended with random numbers previously to make the headers unique.
FormIA_uni_uni <- FormIA_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIA_uni_uni) <- c("Entry",	"Entry name",	"Gene names",	"Gene names  (ORF )",	"Organism",	"Organism ID",	"Length",	"Mass",	"Sequence",	"EC number",	"Fragment",	"Pathway",	"Annotation",	"Status",	"Intramembrane",	"Transmembrane",	"Cross reference (INTERPRO)",	"Cross reference (PFAM)",	"Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)",	"Taxonomic lineage (PHYLUM)",	"Taxonomic lineage (CLASS)",	"Taxonomic lineage (ORDER)",	"Cofactor",	"Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormIA_uni_uni <- FormIA_uni_uni[!is.na(FormIA_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIA_uni_uni) #4 25
head(FormIA_uni_uni)
write.csv(FormIA_uni_uni, "FormIA.csv", row.names=FALSE) 


dim(FormIAc)  #111

FormIAc_uni <- merge(proks_compl, FormIAc, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.

dim(FormIAc_uni)   #179, 26
FormIAc_uni_uni <-  FormIAc_uni[!duplicated(FormIAc_uni[,'Entry']),]
dim(FormIAc_uni_uni)  #49  
FormIAc_uni_uni <- FormIAc_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIAc_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormIAc_uni_uni <- FormIAc_uni_uni[!is.na(FormIAc_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIAc_uni_uni)  #48
head(FormIAc_uni_uni)  #verify that column headers and contents are aligned and properly formatted

write.csv(FormIAc_uni_uni, "FormIAc.csv", row.names=FALSE)



#////////////////////////////////////////////////////////////

# do the rest:
dim(FormIAq) #69 2

FormIAq_uni <- merge(proks_compl, FormIAq, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.
dim(FormIAq_uni)  #113  26
FormIAq_uni_uni <-  FormIAq_uni[!duplicated(FormIAq_uni[,'Entry']),]
dim(FormIAq_uni_uni)  #28 26
FormIAq_uni_uni <- FormIAq_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIAq_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							
FormIAq_uni_uni <- FormIAq_uni_uni[!is.na(FormIAq_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIAq_uni_uni)  #27
head(FormIAq_uni_uni)
write.csv(FormIAq_uni_uni, "FormIAq.csv", row.names=FALSE)


dim(FormIB)  #eukaryotes
dim(proks_compl)
FormIb_uni <- merge(proks_compl, FormIB, by.x = "Sequence", by.y="seq.text", no.dups=FALSE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.
dim(FormIb_uni)
FormIb_uni_uni <-  FormIb_uni[!duplicated(FormIb_uni[,'Entry']),]
dim(FormIb_uni_uni)
#FormIb_uni_uni <- FormIb_uni_uni[ ,c(2:10,1,11:24)]
#colnames(FormIb_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Protein names", "Length", "Mass", "Sequence", "EC number", "Pathway", "Cofactor", "Annotation", "Status", "Subcellular location (CC)", "Intramembrane", "Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)", "Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (KINGDOM)", "Taxonomic lineage (SUBKINGDOM)")
FormIb_uni_uni <- FormIb_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIb_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

head(FormIb_uni_uni)
dim(FormIb_uni_uni)
write.csv(FormIb_uni_uni, "FormIB.csv", row.names=FALSE)

dim(FormIBc)  #103

FormIBc_uni <- merge(proks_compl, FormIBc, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.
dim(FormIBc_uni)  #138 26
FormIBc_uni_uni <-  FormIBc_uni[!duplicated(FormIBc_uni[,'Entry']),]
dim(FormIBc_uni_uni)  #65  26
FormIBc_uni_uni <- FormIBc_uni_uni[ ,c(2:9, 1, 10:25)]
names(FormIBc_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormIBc_uni_uni <- FormIBc_uni_uni[!is.na(FormIBc_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIBc_uni_uni)  #64
head(FormIBc_uni_uni)
write.csv(FormIBc_uni_uni, "FormIBc.csv", row.names=FALSE)

#//////////////////////// troubleshoot unique(), duplicate() for removing repeated protein seq from merged dataframe

#dat <- data.frame(id=c(1,1,3),id2=c(1,1,4),somevalue=c("x","y","z"))
#dat
#dat[!duplicated(dat[,'id']),]  
#dim(dat)  

dim(FormIC)  #112 2
FormIC_uni <- merge(proks_compl, FormIC, by.x = "Sequence", by.y="seq.text", all.y=TRUE) 
dim(FormIC_uni)  #167 26
FormIC_uni_uni <-  FormIC_uni[!duplicated(FormIC_uni[,'Entry']),]
dim(FormIC_uni_uni)  #89 26

FormIC_uni_uni <- FormIC_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIC_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway", "Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormIC_uni_uni <- FormIC_uni_uni[!is.na(FormIC_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIC_uni_uni)  #88 25
head(FormIC_uni_uni)
write.csv(FormIC_uni_uni, "FormIC.csv", row.names=FALSE)

dim(II_III)  #21  2
FormII_III_uni <- merge(proks_compl, II_III, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.
dim(FormII_III_uni)  #21 26

FormII_III_uni_uni <-  FormII_III_uni[!duplicated(FormII_III_uni[,'Entry']),]
dim(FormII_III_uni_uni)  #8 26
FormII_III_uni_uni <- FormII_III_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormII_III_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormII_III_uni_uni <- FormII_III_uni_uni[!is.na(FormII_III_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormII_III_uni_uni)   #7
head(FormII_III_uni_uni)
write.csv(FormII_III_uni_uni, "FormII_III.csv", row.names=FALSE)

dim(FormIE)  #166
FormIe_uni <- merge(proks_compl, FormIE, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all proks_compl rows whose sequences are not found in FormIAc.
dim(FormIe_uni)  #17,000!
#write.csv(FormIe_uni, "FormIe_uni.csv")  #fliterally over 10,000 rows of entry A0A0T9KU63, myucobacterium tuberculosis

FormIe_uni_uni <-  FormIe_uni[!duplicated(FormIe_uni[,'Entry']),]
dim(FormIe_uni_uni) #17  many, many, many mycobacterium tuberculosis entries are possibly duplicates that are removed?
FormIe_uni_uni <- FormIe_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIe_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

#Entry	Entry.name	Gene.names	Gene.names...ORF..	Organism	Organism.ID	Mass	x	Sequence	EC.number	Fragment	Pathway	Cofactor	Annotation	Status	Intramembrane	Transmembrane	Cross.reference..INTERPRO.	Cross.reference..PFAM.	Cross.reference..PROSITE.	Taxonomic.lineage..SUPERKINGDOM.	Taxonomic.lineage..PHYLUM.	Taxonomic.lineage..CLASS.	Taxonomic.lineage..ORDER.	Length	Protein.names																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																						
FormIe_uni_uni <- FormIe_uni_uni[!is.na(FormIe_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIe_uni_uni)  #16
head(FormIe_uni_uni)
write.csv(FormIe_uni_uni, "FormIe.csv", row.names=FALSE)

dim(FormII)  #126
FormII_uni <- merge(proks_compl, FormII, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.
dim(FormII_uni)  #126
FormII_uni_uni <-  FormII_uni[!duplicated(FormII_uni[,'Entry']),]
dim(FormII_uni_uni)  #38
FormII_uni_uni <- FormII_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormII_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormII_uni_uni <- FormII_uni_uni[!is.na(FormII_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormII_uni_uni)  #37
head(FormII_uni_uni)
write.csv(FormII_uni_uni, "FormII.csv", row.names=FALSE)

dim(FormIII)  #40
FormIII_uni <- merge(proks_compl, FormIII, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.
dim(FormIII_uni)  #57
FormIII_uni_uni <-  FormIII_uni[!duplicated(FormIII_uni[,'Entry']),]
dim(FormIII_uni_uni)  #14
FormIII_uni_uni <- FormIII_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIII_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormIII_uni_uni <- FormIII_uni_uni[!is.na(FormIII_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIII_uni_uni)  #36
head(FormIII_uni_uni)
write.csv(FormIII_uni_uni, "FormIII.csv", row.names=FALSE)


dim(FormIIIb)  #97
FormIIIb_uni <- merge(proks_compl, FormIIIb, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.
dim(FormIIIb_uni)  #110

FormIIIb_uni_uni <-  FormIIIb_uni[!duplicated(FormIIIb_uni[,'Entry']),]

dim(FormIIIb_uni_uni)  #82
FormIIIb_uni_uni <- FormIIIb_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIIIb_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway", "Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormIIIb_uni_uni <- FormIIIb_uni_uni[!is.na(FormIIIb_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIIIb_uni_uni)  #81
head(FormIIIb_uni_uni)
write.csv(FormIIIb_uni_uni, "FormIIIb.csv", row.names=FALSE)

dim(FormIIIc)  #2
FormIIIc_uni <- merge(proks_compl, FormIIIc, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.
dim(FormIII_uni)  #57
FormIIIc_uni_uni <-  FormIIIc_uni[!duplicated(FormIIIc_uni[,'Entry']),]
dim(FormIIIc_uni_uni)  #1
head(FormIIIc_uni_uni)
#FormIII_uni_uni <- FormIII_uni_uni[ ,c(2:10,1,11:24)]
#colnames(FormIII_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Protein names", "Length", "Mass", "Sequence", "EC number", "Pathway", "Cofactor", "Annotation", "Status", "Subcellular location (CC)", "Intramembrane", "Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)", "Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (KINGDOM)", "Taxonomic lineage (SUBKINGDOM)")
FormIIIc_uni_uni <- FormIIIc_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIIIc_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormIIIc_uni_uni <- FormIIIc_uni_uni[!is.na(FormIIIc_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIIIc_uni_uni)

head(FormIIIc_uni_uni)

names(FormIIIc_uni_uni)
write.csv(FormIIIc_uni_uni, "FormIIIc.csv", row.names=FALSE)

dim(FormIV)  #1076
FormIV_uni <- merge(proks_compl, FormIV, by.x = "Sequence", by.y="seq.text", all.y=TRUE)  #this should remove all rub_uniprot rows whose sequences are not found in FormIAc.

dim(FormIV_uni)  #20,054
FormIV_uni_uni <-  FormIV_uni[!duplicated(FormIV_uni[,'Entry']),]
dim(FormIV_uni_uni)  #457
FormIV_uni_uni <- FormIV_uni_uni[ ,c(2:9,1,10:25)]  #re-arrange columns to put them in same order as uniprot downloads. remove the fragment column because its not found in the other .csvs at this point.
names(FormIV_uni_uni) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "Sequence", "EC number", "Fragment", "Pathway",	"Annotation", "Status", "Intramembrane",	"Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)",	"Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Cofactor", "Protein names")																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							

FormIV_uni_uni <- FormIV_uni_uni[!is.na(FormIV_uni_uni$Entry),]  #get rid of all rows in df that have NA in the "Entry" column. 
dim(FormIV_uni_uni)  #456

head(FormIV_uni_uni)

write.csv(FormIV_uni_uni, "FormIV.csv", row.names=FALSE)

#//////////////////////////////////////////////

#coooool ok I want to make more files, FormIall, FormIIIall, etc.



FormI_all <- rbind(FormI_uni_uni, FormIA_uni_uni, FormIAc_uni_uni, FormIAq_uni_uni, FormIBc_uni_uni, FormIC_uni_uni, FormIe_uni_uni)
head(FormI_uni_uni)
head(FormIA_uni_uni)
dim(FormI_all) #225
head(FormIAc_uni_uni)
#FormI_all <- rbind(FormI_uni_uni, FormIA_uni_uni, FormIAc_uni_uni, FormIAq_uni_uni, FormIBc_uni_uni, FormIC_uni_uni, FormIe_uni_uni)  #no form Id; these are euks.
dim(FormI_all)
names(FormI_all)
head(FormI_all)
#colnames(FormI_all) <- c("Entry", "Entry name", "Gene names", "Gene names  (ORF )", "Organism", "Organism ID", "Length", "Mass", "X", "Sequence", "EC number", "Pathway", "Cofactor", "Annotation", "Status", "Intramembrane", "Transmembrane", "Cross reference (INTERPRO)", "Cross reference (PFAM)", "Cross reference (PROSITE)", "Taxonomic lineage (SUPERKINGDOM)", "Taxonomic lineage (PHYLUM)", "Taxonomic lineage (CLASS)", "Taxonomic lineage (ORDER)", "Protein names")

names(FormI_all)
write.csv(FormI_all, "FormI_all.csv", row.names=FALSE)

FormIA_all <- rbind(FormIA_uni_uni, FormIAc_uni_uni, FormIAq_uni_uni)
head(FormIA_all)
write.csv(FormIA_all, "FormIA_all.csv", row.names=FALSE)


FormII_all <- rbind(FormII_III_uni_uni, FormII_uni_uni)
head(FormII_all)
dim(FormII_all)
write.csv(FormII_all, "FormII_all.csv", row.names=FALSE)

formIII_all <- rbind(FormIII_uni_uni, FormIIIb_uni_uni, FormIIIc_uni_uni)
dim(formIII_all)
head(formIII_all)
write.csv(formIII_all, "FormIII_all.csv", row.names=FALSE)

dim(formIII_all)
dim(FormI_uni_uni)
dim(FormIAc_uni_uni)
dim(FormIAq_uni_uni)
dim(FormIBc_uni_uni)
dim(FormIb_uni_uni)
dim(FormIC_uni_uni)
dim(FormIe_uni_uni)
dim(FormII_uni_uni)
dim(FormIII_uni_uni)
dim(FormIV_uni_uni)
dim(FormI_all)
getwd()
#names(FormI_all)
#names(FormII_uni_uni)
#dim(FormIAc_uni_uni)

#Rub_lg_all <- rbind(FormI_all, FormII_uni_uni, FormIII_uni_uni, FormIV_uni_uni)
#dim(FormI_all)

###### now create files for oxyphen input:

setwd("/Users/tonglen/Downloads/oxyphen-multiome/Org_ID_inputs")

FormIBc_org_id <- FormIBc_uni_uni[ ,c(6, 1:5, 7:25)]
names(FormIBc_org_id)[1] <- "Organism_ID"
head(FormIBc_org_id)
write.csv(FormIBc_org_id, "FormIBc_org_id.csv", row.names=FALSE)

FormII_all_org_id <- FormII_all[ ,c(6, 1:5, 7:25)]
names(FormII_all_org_id)[1] <- "Organism_ID"
head(FormII_all_org_id)
write.csv(FormII_all_org_id, "FormII_all_org_id.csv", row.names=FALSE)

FormIIIb_org_id <- FormIIIb_uni_uni[ ,c(6, 1:5, 7:25)]
names(FormIIIb_org_id)[1] <- "Organism_ID"
head(FormIIIb_org_id)
write.csv(FormIIIb_org_id, "FormIIIb_org_id.csv", row.names=FALSE)



#/////////////////////////

#avi's code is choking on the new files, and the only thing
# i think it could be is the order of columns is different in my files.
# so, try re-ordering the columns in one file and check.
# the order should be:
#names(rub_uniprot)
#names(FormI_all)

#names(FormI_all)

#/////////////////////////////////

# troubleshooting some bullshit with avi's code re: jaccard index

#setwd("/Users/tonglen/avi script/rub-co-occurrence-master/UniprotKB_downloads/")

#bin_proks <- read.table("binary_prokaryote_df.csv", sep = '\t', header = T)

#intersection <- as.data.frame(intersect(bin_proks$Organism_ID, FormI_all$Organism.ID))

#dim(intersection)
#head(intersection)

#head(FormI_uni_uni$Sequence, 100)
