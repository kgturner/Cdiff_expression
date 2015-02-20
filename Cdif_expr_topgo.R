#C.diffusa expression - topGO
#2/17/15
#with help from Kay, see email chain "questions about microarray analysis" starting 11/14/14

source("http://bioconductor.org/biocLite.R")
biocLite(c("topGO", "ALL", "genefilter", "multtest"))

library(topGO)
# library("ALL")

#load data
PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") #contigs and q values
PC1p <- read.table("lme4_PC1_pvalues.txt", header=T, sep="\t") #contigs and p values
PC1q_intsigdf <- read.table("PC1_sigint_df.txt", header=T)

#make custom annotation obj
# TR001 <- read.table("~/GitHub/Cdiff_expression/GOanalysis/unique.sorted.out.ath_Cendif1.unigenes_GO_flip", header=T) 
# US022 <- read.table("~/GOanalysis/unique.sorted.out.ath_Cendif2.unigenes_GO_flip", header=T)

# TR001 <- readMappings(file = "~/GitHub/Cdiff_expression/GOanalysis/unique.sorted.out.ath_Cendif1.unigenes_GO_flip") #go annots from BLAST
# US022 <- readMappings(file = "~/GitHub/Cdiff_expression/GOanalysis/unique.sorted.out.ath_Cendif2.unigenes_GO_flip")
GOmap <- readMappings(file = "~/GitHub/Cdiff_expression/GOanalysis/GOmap_Cendif1_Cendif2.txt_flip")
str(head(GOmap))

#all contigs, i.e. gene universe
CdifNames <- names(GOmap)
head(CdifNames)

#identify genes of interest
# #list of sig/not siq
# intqList <- factor(as.integer(PC1q$intQsig))
# names(intqList) <- CdifNames
# str(intqList)

#OR, incorporating q-values for additional analysis options
intqList <- PC1q$intQ
names(intqList) <- rownames(PC1q)
#set alpha level
topDiffGenes <- function(allScore) {
  return(allScore < 0.05)
}

#make topGOdata object
intqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                  ontology="BP", #MF, CC?
                  allGenes=intqList, #factor describing which genes are of interest/sig, which are not
                  geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                  annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                  nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                  gene2GO=GOmap) #our gene->GO term mapping file

intqGOdata
#stats for random sample
sel.terms <- sample(usedGO(intqGOdata), 10)
termStat(intqGOdata, sel.terms)

#enrichment
#kay used fisher exact test...
resultFis <- runTest(intqGOdata, algorithm = "classic", statistic = "fisher")

pvalFis <- score(resultFis)
head(pvalFis)
hist(pvalFis, 50, xlab = "p-values")

allRes <- GenTable(intqGOdata, classic = resultFis, orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj


####example####
data(ALL) #eset
data(geneList) #go term and p val
affyLib <- paste(annotation(ALL), "db", sep = ".")
biocLite("hgu95av2.db")
library(package = affyLib, character.only = TRUE)
sum(topDiffGenes(geneList))

geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
head(geneID2GO)

geneNames <- names(geneID2GO)
head(geneNames)
myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))

names(geneList) <- geneNames
str(geneList)

#version w/ pvals?
library(genefilter)
library(multtest)
selProbes <- genefilter(ALL, filterfun(pOverA(0.20, log2(100)), function(x) (IQR(x) > 0.25)))
eset <- ALL[selProbes, ]
y <- as.integer(sapply(eset$BT, function(x) return(substr(x, 1, 1) == 'T')))
table(y)
geneList <- getPvalues(exprs(eset), classlabel = y, alternative = "greater")


topDiffGenes <- function(allScore) {
  return(allScore < 0.05)
  }
x <- topDiffGenes(geneList)
sum(x) ## the number of selected genes
GOdata <- new("topGOdata",
              + description = "GO analysis of ALL data; B-cell vs T-cell",
              + ontology = "BP",
              + allGenes = geneList,
              + geneSel = topDiffGenes,
              + annot = annFUN.db,
              + nodeSize = 5,
              + affyLib = affyLib)