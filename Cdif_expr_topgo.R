#C.diffusa expression - topGO
#2/17/15
#from Kay, see email chain "questions about microarray analysis" starting 11/14/14

source("http://bioconductor.org/biocLite.R")
biocLite(c("topGO", "ALL"))

library(topGO)
# library("ALL")

#load data
PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") #contigs and p and q values
PC1q_intsigdf <- read.table("PC1_sigint_df.txt", header=T)
TR001 <- read.table("~/GOanalysis/unique.sorted.out.ath_Cendif1.unigenes_GO_flip", header=T) #go annots from BLAST
US022 <- read.table("~/GOanalysis/unique.sorted.out.ath_Cendif2.unigenes_GO_flip", header=T)

#subset PC1q to sig int, merge w/ go annoations
intq <- subset(PC1q,intQsig==TRUE, select=c("Contig", "intQ") ) #contigs with sig Origin*Trt

#example
data(ALL) #eset
data(geneList) #go term and p val
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)
