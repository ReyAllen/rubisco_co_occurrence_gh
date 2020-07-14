# purpose of this code is to make pretty cdf of oxyphen data in ggplot2, hopefully,
# similar to the figs in the original oxyphen paper (stepped line graphs).

#step one: import dataframe of stuff
#gonna want plots at 60% identity and at 40.  Prob start with 60. 
#step two: define factors for plottinggggg
#step three: plot shit.

#Re-used loop from oxyphen code write fasta from csv

# hadley wickam being a dick on a public forum: https://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-data-frames
rm(list=ls())
setwd("~/Downloads/oxyphen-multiome/OUTPUT_60/")  #this folder contains multiple folders with the 60% homology oxyphen data, one Rubisco isoform per _60 folder.

library(phylotools)
library(tidyr)
library(dplyr)

filenames <- list.files(path="/Users/tonglen/Downloads/oxyphen-multiome/OUTPUT_60/", pattern="*.tsv", full.names=FALSE, recursive=FALSE)
head(filenames)

test <- read.csv("results_table_IBc_60.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
head(test)

##Create list of data frame names without the ".csv" part 

df_names <- list()    #make an empty list to store names of new dataframes
for (i in 1:length(filenames))  #iterate through the filenames list
{  
  handle <- unlist(strsplit(filenames[i], "\\."))  #split filenames into filename and .tsv
  handle <- handle[-2]                              #remove .tsv from vector filenames[i]
  handle <-substring(handle,15)
  df_names <- append(df_names, handle)                  #put your most recent new filename in the list of names you made
}
df_names


for(i in df_names)
  {
  filepath <- file.path("/Users/tonglen/Downloads/oxyphen-multiome/OUTPUT_60/", paste("results_table_", i, ".tsv",sep=""))
  assign (i, read.delim(filepath, colClasses=c("character","numeric","character"), sep = "\t", header = FALSE))
  }

head(IBc_60)  #woo there it iiiiiiiiis >
# ////////////////////  Makin' the graphs /////////////////////////

#http://www.sthda.com/english/wiki/ggplot2-ecdf-plot-quick-start-guide-for-empirical-cumulative-density-function-r-software-and-data-visualization

library(ggplot2)

#make dataframe for plotting using ggplot2



Isoform <- c(rep("Ru. IBc", length(IBc_60$V2)),  #tag is for creating legend in ggplot2
           rep("Ru. II", length(II_all_60$V2)),
           rep("Ru. IIIb", length(IIIb_60$V2)))

df <- data.frame(c(IBc_60$V2, II_all_60$V2, IIIb_60$V2), Isoform) #create df of two columns: first column data, second column data label
dim(df)
head(df)

colnames(df) <- c("Counts", "Isoform")
head(df)

#remove data from 1904441, which has 48 o2-utilzying enzymes (throws off x-axis), and is uncharacterized alphaproteob.
#https://www.uniprot.org/taxonomy/1904441

dim(df) #136
df <- df[df$Counts != 48, ]
dim(df) #135

sixty <- ggplot(df, aes(df$Counts, color = Isoform)) + stat_ecdf(geom = "step")+
  labs(title="60% identity to O2-utilizing enzymes",
       y = "cumulative density", x="No. of O2-utilizing enzyme classes per proteome")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
  
sixty  

#add legend with entries from multiple sources:
#https://aosmith.rbind.io/2018/07/19/manual-legends-ggplot2/
#geom_smooth(method = "lm", se = FALSE, aes(color = "black") )


sixty + geom_vline(xintercept=2, linetype="dashed", size=0.5, color="grey")+
  geom_vline(xintercept=9, linetype="dotted", size=0.5, color="grey")
#  scale_color_identity(guide = "legend")
#geom_vline(aes(xintercept=years, color=labels),data=vertDf, show_guide=T)

#////////////////////////// plot below using base functions used to work but doesn't anymore; legend takes over whole plot.

# https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics

par(pty='s')  # force the plot to be square before we start
#par(xpd=FALSE)  #turn on clipping

plot(ecdf(IBc_60$V2),
     xlim=c(0,30),
     xlab="No. of O2-utilizing enzymes per proteome",
     ylab="Cumulative density",
     main="60% identity to O2-utilizing enzymes",
     col="cornflowerblue",
     pch = 15, 
     frame.plot=FALSE,
     cex=1)

lines(ecdf(IIIb_60$V2), pch = 15, cex = 1,
      col="coral")
lines(ecdf(II_all_60$V2), pch = 15, cex = 1,
      col="chartreuse")
abline(v=9, lty = "dotted", #make into vertical dashed line for aerobe cutoff
      col="deepskyblue")
abline(v=2, lty = "dashed",   #make into vertical dashed line for anaerobe cutoff
       col = "orange")
#par(xpd=TRUE) #turn off clipping before making legend
legend('bottomright',
      legend=c("Ru.IBc","Ru.IIIb","Ru.II","Aerobe", "Anaerobe"), #text in legend
      col=c("cornflowerblue","coral","chartreuse","deepskyblue", "orange"),  # point colors
      pch=15, # specify the point type to be a square
      cex = 0.4 #relative size of legend
      )  

dim(IBc_60)  #get sample sizes, n approx. 200 per isoform for 60% identity analyses
dim(IIIb_60)
dim(II_all_60)

#////////////////  repeat all of above for 40% identity data //////////////////////////////////////


setwd("~/Downloads/oxyphen-multiome/OUTPUT_40/")  #this folder contains multiple folders with the 60% homology oxyphen data, one Rubisco isoform per _60 folder.

library(phylotools)
library(tidyr)
library(dplyr)

filenames <- list.files(path="/Users/tonglen/Downloads/oxyphen-multiome/OUTPUT_40/", pattern="*.tsv", full.names=FALSE, recursive=FALSE)
head(filenames)

test <- read.csv("results_table_IBc_40.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
head(test)

##Create list of data frame names without the ".tsv" part 

df_names <- list()    #make an empty list to store names of new dataframes
for (i in 1:length(filenames))  #iterate through the filenames list
{  
  handle <- unlist(strsplit(filenames[i], "\\."))  #split filenames into filename and .tsv
  handle <- handle[-2]                              #remove .tsv from vector filenames[i]
  handle <-substring(handle,15)
  df_names <- append(df_names, handle)                  #put your most recent new filename in the list of names you made
}
df_names


for(i in df_names)
{
  filepath <- file.path("/Users/tonglen/Downloads/oxyphen-multiome/OUTPUT_40/", paste("results_table_", i, ".tsv",sep=""))
  assign (i, read.delim(filepath, colClasses=c("character","numeric","character"), sep = "\t", header = FALSE))
}

head(IBc_40)  #woo there it iiiiiiiiis E>

# ////////////////////  Makin' the graphs /////////////////////////

head(IBc_40$V2)


par(pty='s')  # force the plot to be square before we start

plot(ecdf(IBc_40$V2),
     xlim=c(0,50),
     xlab="No. of O2-utilizing EC classes per proteome",
     ylab="Cumulative density",
     main="40% identity to O2-utilizing enzymes",
     col="cornflowerblue",
     pch = 15, 
     frame.plot=FALSE,
     cex=1)

lines(ecdf(IIIa_40$V2), pch = 15, cex = 1,
      col="coral")
lines(ecdf(II_40$V2), pch = 15, cex = 1,
      col="chartreuse")
abline(v=29, lty = "dotted", #make into vertical dashed line for aerobe cutoff
       col="deepskyblue")
abline(v=13, lty = "dashed",   #make into vertical dashed line for anaerobe cutoff
       col = "orange")
legend('bottomright', 
       legend=c("Ru.IB","Ru.IIIa","Ru.II","Aerobe", "Anaerobe"),  # text in the legend
       col=c("cornflowerblue","coral","chartreuse","deepskyblue", "orange"),  # point colors
       pch=15)  # specify the point type to be a square

dim(IBc_40)
dim(II_40)
dim(IIIa_40)
