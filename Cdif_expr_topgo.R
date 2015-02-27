#C.diffusa expression - topGO
#2/17/15
#with help from Kay, see email chain "questions about microarray analysis" starting 11/14/14

source("http://bioconductor.org/biocLite.R")
biocLite(c("topGO", "ALL", "genefilter", "multtest"))

library(topGO)
# library("ALL")
library(qvalue)

#load data
PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") #contigs and q values
PC1p <- read.table("lme4_PC1_pvalues.txt", header=T, sep="\t") #contigs and p values
PC1q_intsigdf <- read.table("PC1_sigint_df.txt", header=T)
PC1q_Osigdf <- read.table("PC1_sigOrigin_df.txt", header=T)
PC1q_sigdf <- read.table("PC1_sigPC1_df.txt", header=T)
PC1q_trtsigdf <- read.table("PC1_sigtrt_df.txt", header=T)

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

#set alpha level
topDiffGenes <- function(allScore) {
  return(allScore < 0.05)
}

####int sig####
#identify genes of interest
# #list of sig/not siq
# intqList <- factor(as.integer(PC1q$intQsig))
# names(intqList) <- CdifNames
# str(intqList)

#OR, incorporating q-values for additional analysis options
intqList <- PC1q$intQ
names(intqList) <- PC1q$Contig
head(intqList)

#make topGOdata object
intqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                  ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                  allGenes=intqList, #factor describing which genes are of interest/sig, which are not
                  geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                  annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                  nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                  gene2GO=GOmap) #our gene->GO term mapping file

intqGOdata
# #stats for random sample
# sel.terms <- sample(usedGO(intqGOdata), 10)
# termStat(intqGOdata, sel.terms)

#enrichment
#kay used fisher exact test...
resultclas <- runTest(intqGOdata, algorithm = "classic", statistic = "fisher")
resultelim <- runTest(intqGOdata, algorithm = "elim", statistic = "fisher")
resultwt <- runTest(intqGOdata, algorithm = "weight", statistic = "fisher")
resultwt01 <- runTest(intqGOdata, algorithm = "weight01", statistic = "fisher")
resultlea <- runTest(intqGOdata, algorithm = "lea", statistic = "fisher")
resultpc <- runTest(intqGOdata, algorithm = "parentchild", statistic = "fisher")

# pvalFis <- score(resultFis)
# head(pvalFis)
# hist(pvalFis, 50, xlab = "p-values")

#work on this....
clasRes <- as.data.frame(GenTable(intqGOdata, classic = resultclas, topNodes=2483))
clasRes <- GenTable(intqGOdata, classic = resultclas, topNodes=50)
#correct for multiple tests?
clasResQ <- qvalue(p=clasRes$classic, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
               smooth.df=3, smooth.log.pi0=FALSE)
clasRes <- cbind(clasRes, GOtermQ=clasResQ$qvalues, GOtermsig=clasResQ$significant) #intQ=intQ$qvalues, intQsig=intQ$significant
write.table(clasRes, file="GOresults_sigint_classic.txt", sep="\t")

#
allRes <- GenTable(intqGOdata, classic = resultclas, elim=resultelim, weight=resultwt, 
                   weight01=resultwt01, lea=resultlea, parentchild=resultpc,
                   orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj
write.table(allRes, file="GOresults_sigint.txt", sep="\t")

####O sig####
#OR, incorporating q-values for additional analysis options
oqList <- PC1q$originQ
names(oqList) <- PC1q$Contig

#make topGOdata object
oqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin effect",
                  ontology="BP", #MF, CC?
                  allGenes=oqList, #factor describing which genes are of interest/sig, which are not
                  geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                  annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                  nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                  gene2GO=GOmap) #our gene->GO term mapping file

oqGOdata
#stats for random sample
sel.terms <- sample(usedGO(oqGOdata), 10)
termStat(oqGOdata, sel.terms)

#enrichment
#kay used fisher exact test...
resultclas <- runTest(oqGOdata, algorithm = "classic", statistic = "fisher")
resultelim <- runTest(oqGOdata, algorithm = "elim", statistic = "fisher")
resultwt <- runTest(oqGOdata, algorithm = "weight", statistic = "fisher")
resultwt01 <- runTest(oqGOdata, algorithm = "weight01", statistic = "fisher")
resultlea <- runTest(oqGOdata, algorithm = "lea", statistic = "fisher")
resultpc <- runTest(oqGOdata, algorithm = "parentchild", statistic = "fisher")

allResO <- GenTable(oqGOdata, classic = resultclas, elim=resultelim, weight=resultwt, 
                   weight01=resultwt01, lea=resultlea, parentchild=resultpc,
                   orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj
write.table(allResO, file="GOresults_sigOrigin.txt", sep="\t")

####trt sig####
#OR, incorporating q-values for additional analysis options
trtqList <- PC1q$trtQ
names(trtqList) <- PC1q$Contig

#make topGOdata object
trtqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Trt effect",
                  ontology="BP", #MF, CC?
                  allGenes=trtqList, #factor describing which genes are of interest/sig, which are not
                  geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                  annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                  nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                  gene2GO=GOmap) #our gene->GO term mapping file

trtqGOdata
# #stats for random sample
# sel.terms <- sample(usedGO(trtqGOdata), 10)
# termStat(trtqGOdata, sel.terms)

#enrichment
#kay used fisher exact test...
resultclas <- runTest(trtqGOdata, algorithm = "classic", statistic = "fisher")
resultelim <- runTest(trtqGOdata, algorithm = "elim", statistic = "fisher")
resultwt <- runTest(trtqGOdata, algorithm = "weight", statistic = "fisher")
resultwt01 <- runTest(trtqGOdata, algorithm = "weight01", statistic = "fisher")
resultlea <- runTest(trtqGOdata, algorithm = "lea", statistic = "fisher")
resultpc <- runTest(trtqGOdata, algorithm = "parentchild", statistic = "fisher")

allResTrt <- GenTable(trtqGOdata, classic = resultclas, elim=resultelim, weight=resultwt, 
                   weight01=resultwt01, lea=resultlea, parentchild=resultpc,
                   orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj
write.table(allResTrt, file="GOresults_sigTrt.txt", sep="\t")

####PC1 sig####
#OR, incorporating q-values for additional analysis options
pcqList <- PC1q$covQ
names(pcqList) <- PC1q$Contig

#make topGOdata object
pcqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig PC1 effect",
                  ontology="BP", #MF, CC?
                  allGenes=pcqList, #factor describing which genes are of interest/sig, which are not
                  geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                  annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                  nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                  gene2GO=GOmap) #our gene->GO term mapping file

pcqGOdata
# #stats for random sample
# sel.terms <- sample(usedGO(intqGOdata), 10)
# termStat(intqGOdata, sel.terms)

#enrichment
#kay used fisher exact test...
resultclas <- runTest(pctqGOdata, algorithm = "classic", statistic = "fisher")
resultelim <- runTest(pcqGOdata, algorithm = "elim", statistic = "fisher")
resultwt <- runTest(pcqGOdata, algorithm = "weight", statistic = "fisher")
resultwt01 <- runTest(pcqGOdata, algorithm = "weight01", statistic = "fisher")
resultlea <- runTest(pcqGOdata, algorithm = "lea", statistic = "fisher")
resultpc <- runTest(pcqGOdata, algorithm = "parentchild", statistic = "fisher")

allResPC1 <- GenTable(pcqGOdata, classic = resultclas, elim=resultelim, weight=resultwt, 
                   weight01=resultwt01, lea=resultlea, parentchild=resultpc,
                   orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj
write.table(allResPC1, file="GOresults_sigPC1.txt", sep="\t")

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