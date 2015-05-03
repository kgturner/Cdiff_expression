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
# PC1p <- read.table("lme4_PC1_pvalues.txt", header=T, sep="\t") #contigs and p values
# PC1q_intsigdf <- read.table("PC1_sigint_df.txt", header=T)
# PC1q_Osigdf <- read.table("PC1_sigOrigin_df.txt", header=T)
# PC1q_sigdf <- read.table("PC1_sigPC1_df.txt", header=T)
# PC1q_trtsigdf <- read.table("PC1_sigtrt_df.txt", header=T)

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

##########################################
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

#make topGOdata object, BP
intqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                  ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                  allGenes=intqList, #factor describing which genes are of interest/sig, which are not
                  geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                  annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                  nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                  gene2GO=GOmap) #our gene->GO term mapping file

intqGOdata

#   
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 61024 available genes (all genes from the array):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  0.7364 0.761989 0.5792765 0.78036662 0.471659746  ...
# - 227  significant genes. 
# 
# 37373 feasible genes (genes that can be used in the analysis):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  0.7364 0.761989 0.5792765 0.78036662 0.471659746  ...
# - 180  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2483 
# - number of edges = 5282 



#enrichment
#kay used fisher exact test...
# resultclas <- runTest(intqGOdata, algorithm = "classic", statistic = "fisher")
# resultelim <- runTest(intqGOdata, algorithm = "elim", statistic = "fisher")
# resultwt <- runTest(intqGOdata, algorithm = "weight", statistic = "fisher")
# resultwt01 <- runTest(intqGOdata, algorithm = "weight01", statistic = "fisher")
# resultlea <- runTest(intqGOdata, algorithm = "lea", statistic = "fisher")

resultpc <- runTest(intqGOdata, algorithm = "parentchild", statistic = "fisher")

# 
# allRes <- GenTable(intqGOdata, classic = resultclas, elim=resultelim, weight=resultwt, 
#                    weight01=resultwt01, lea=resultlea, parentchild=resultpc,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj
# write.table(allRes, file="GOresults_sigint.txt", sep="\t")

#multiple testing correction needed?
# pvalFis <- score(resultclas)
# head(pvalFis)
# hist(pvalFis, 50, xlab = "p-values")
# #correct for multiple tests?
# clasResQ <- qvalue(p=pvalFis, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
#                smooth.df=3, smooth.log.pi0=FALSE)
# clasRes <- cbind(clasRes, GOtermQ=clasResQ$qvalues, GOtermsig=clasResQ$significant) #intQ=intQ$qvalues, intQsig=intQ$significant
# write.table(clasRes, file="GOresults_sigint_classic.txt", sep="\t")
# 
# or from topGO documentation:"For the methods that account for the GO topology like elim and weight...multiple testing theory does not
# directly apply. We like to interpret the p-values returned by these methods as corrected or not affected
# by multiple testing." 
# topoRes <- GenTable(intqGOdata, elim = resultelim, weight = resultwt, orderby="weight", ranksof="elim", topNodes=50)

#explaining GenTable output, from https://stat.ethz.ch/pipermail/bioconductor/2009-July/028616.html:
# The "Annotated", "Significant" and "Expected" columns show
# statistics computed for each GO term based on the complete
# annotations, meaning that the "true path rule" is used to annotated
# the genes to higher level terms. The "Expected" column shows an
# estimate of the number of genes, anode of size "Annotated" will have
# if the significant genes would be randomly selected from the gene
# universe. Now, if you would use the "classic" algorithm for testing
# for over-representation, then all GO terms with significant values
# will have the "Significant"  < "Expected". However this is not the
# case when using methods like "elim" or "weight" which remove or weight
# genes annotated to GO terms when computing the significance. This
# happens because when you "remove" the genes the counts for the
# specific GO term change and the ratio between "Significant" and
# "Expected" changes. 

# Used parentchild in InvSyn paper, so...
resultpc
# 
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2483 GO terms scored: 33 terms with p < 0.01
# Annotation data:
#   Annotated genes: 37373 
# Significant genes: 180 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 639 

pcRes <- GenTable(intqGOdata, parentchild = resultpc, topNodes=33)
write.table(pcRes, file="GOresults_sigint_pc.txt", sep="\t")

# get more info on specific terms, here, top 3 
mget(pcRes[1:3,1], GOTERM)


####int sig, up in inv, tmpt 2, drought####
int2dsummary <- read.table(file="sigint_popMeans_sumT2drought.txt", header=T)
summary(int2dsummary)
invUp2dList <- subset(int2dsummary, InvUp==TRUE)
invUp2dList.1 <- invUp2dList$Contig

#identify genes of interest
# #list of sig/not siq
int_invUp2dList <- factor(as.integer(CdifNames %in% invUp2dList.1))
names(int_invUp2dList) <- CdifNames
str(int_invUp2dList)
summary(int_invUp2dList)

invUp2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                  ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                  allGenes=int_invUp2dList, #factor describing which genes are of interest/sig, which are not
#                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                  annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                  nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                  gene2GO=GOmap) #our gene->GO term mapping file
invUp2dGOdata
------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 146  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 137  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 
# 
------------------------- topGOdata object -------------------------
  
resultpc <- runTest(invUp2dGOdata, algorithm = "parentchild", statistic = "fisher")
resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 32 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 137 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 516 

invUp2d_pcRes <- GenTable(invUp2dGOdata, parentchild = resultpc, topNodes=32)
write.table(invUp2d_pcRes, file="GOresults_sigint_invUp2drought_pc.txt", sep="\t")

invUp2d_pcRes.more <- GenTable(invUp2dGOdata, parentchild = resultpc, topNodes=70)

# get more info on specific terms, here, top 3 
mget(pcRes[1:3,1], GOTERM)



####int sig, down in inv, tmpt 2, drought####
int2dsummary <- read.table(file="sigint_popMeans_sumT2drought.txt", header=T)

invDn2dList <- subset(int2dsummary, InvUp==FALSE)
invDn2dList.1 <- invDn2dList$Contig

#identify genes of interest
# #list of sig/not siq
int_invDn2dList <- factor(as.integer(CdifNames %in% invDn2dList.1))
names(int_invDn2dList) <- CdifNames
str(int_invDn2dList)
summary(int_invDn2dList)

invDn2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_invDn2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
invDn2dGOdata
------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 49  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 43  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 
# 
------------------------- topGOdata object -------------------------
#   
  
resultpc <- runTest(invDn2dGOdata, algorithm = "parentchild", statistic = "fisher")
resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 8 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 43 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 289 

invDn2d_pcRes <- GenTable(invDn2dGOdata, parentchild = resultpc, topNodes=8)
write.table(invDn2d_pcRes, file="GOresults_sigint_invDn2drought_pc.txt", sep="\t")

invDn2d_pcRes.more<- GenTable(invDn2dGOdata, parentchild = resultpc, topNodes=50)

# get more info on specific terms, here, top 3 
# mget(pcRes[1:3,1], GOTERM)


####int sig, up in inv, tmpt 2 control####
int2csummary <- read.table(file="sigint_popMeans_sumT2Control.txt", header=T)

invUp2cList <- subset(int2csummary, InvUp==TRUE)
invUp2cList.1 <- invUp2cList$Contig

#identify genes of interest
# #list of sig/not siq
int_invUp2cList <- factor(as.integer(CdifNames %in% invUp2cList.1))
names(int_invUp2cList) <- CdifNames
str(int_invUp2cList)
summary(int_invUp2cList)

invUp2cGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_invUp2cList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
invUp2cGOdata
------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 60  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 51  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 

------------------------- topGOdata object -------------------------

resultpc <- runTest(invUp2cGOdata, algorithm = "parentchild", statistic = "fisher")
resultpc

# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 6 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 51 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 371 

invUp2c_pcRes <- GenTable(invUp2cGOdata, parentchild = resultpc, topNodes=6)
write.table(invUp2c_pcRes, file="GOresults_sigint_invUp2control_pc.txt", sep="\t")

invUp2c_pcRes.more <- GenTable(invUp2cGOdata, parentchild = resultpc, topNodes=30)


####int sig, down in inv, tmpt2 control####
int2csummary <- read.table(file="sigint_popMeans_sumT2Control.txt", header=T)

invDn2cList <- subset(int2csummary, InvUp==FALSE)
invDn2cList.1 <- invDn2cList$Contig

#identify genes of interest
# #list of sig/not siq
int_invDn2cList <- factor(as.integer(CdifNames %in% invDn2cList.1))
names(int_invDn2cList) <- CdifNames
str(int_invDn2cList)
summary(int_invDn2cList)

invDn2cGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_invDn2cList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
invDn2cGOdata
------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 135  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 129  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 

------------------------- topGOdata object -------------------------
  
resultpc <- runTest(invDn2cGOdata, algorithm = "parentchild", statistic = "fisher")
resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 34 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 129 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 481 

invDn2c_pcRes <- GenTable(invDn2cGOdata, parentchild = resultpc, topNodes=34)
write.table(invDn2c_pcRes, file="GOresults_sigint_invDn2control_pc.txt", sep="\t")

invDn2c_pcRes.more <- GenTable(invDn2cGOdata, parentchild = resultpc, topNodes=75)

####MF int sig####
# intqList <- PC1q$intQ
# names(intqList) <- PC1q$Contig
# head(intqList)
#make topGOdata object, MF
MFintqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                    ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                    allGenes=intqList, #factor describing which genes are of interest/sig, which are not
                    geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                    annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                    nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                    gene2GO=GOmap) #our gene->GO term mapping file

MFintqGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 61024 available genes (all genes from the array):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  0.7364 0.761989 0.5792765 0.78036662 0.471659746  ...
# - 227  significant genes. 
# 
# 37190 feasible genes (genes that can be used in the analysis):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  0.7364 0.761989 0.5792765 0.78036662 0.471659746  ...
# - 187  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1100 
# - number of edges = 1376 

#enrichment
MFresultpc <- runTest(MFintqGOdata, algorithm = "parentchild", statistic = "fisher")
MFresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1100 GO terms scored: 7 terms with p < 0.01
# Annotation data:
#   Annotated genes: 37190 
# Significant genes: 187 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 239 

(MFpcRes <- GenTable(MFintqGOdata, parentchild = MFresultpc, topNodes=7))
# MFpcRes
write.table(MFpcRes, file="GOresults_MFsigint_pc.txt", sep="\t")

(MFpcRes.more <- GenTable(MFintqGOdata, parentchild = MFresultpc, topNodes=23))

####MF int sig, up in inv, tmpt 2, drought####
int2dsummary <- read.table(file="sigint_popMeans_sumT2drought.txt", header=T)
# summary(int2dsummary)
invUp2dList <- subset(int2dsummary, InvUp==TRUE)
invUp2dList.1 <- invUp2dList$Contig

#identify genes of interest
# #list of sig/not siq
int_invUp2dList <- factor(as.integer(CdifNames %in% invUp2dList.1))
names(int_invUp2dList) <- CdifNames
# str(int_invUp2dList)
# summary(int_invUp2dList)

MFinvUp2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=int_invUp2dList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
MFinvUp2dGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 146  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 139  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635 
MFinvUp2d_resultpc <- runTest(MFinvUp2dGOdata, algorithm = "parentchild", statistic = "fisher")
MFinvUp2d_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 7 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 139 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 194 
MFinvUp2d_pcRes <- GenTable(MFinvUp2dGOdata, parentchild = MFinvUp2d_resultpc, topNodes=7)
MFinvUp2d_pcRes
write.table(MFinvUp2d_pcRes, file="GOresults_MFsigint_invUp2drought_pc.txt", sep="\t")
# 
# invUp2d_pcRes.more <- GenTable(invUp2dGOdata, parentchild = resultpc, topNodes=70)
# 
mget(MFinvUp2d_pcRes[,1], GOTERM)

####MF int sig, down in inv, tmpt 2, drought####
# int2dsummary <- read.table(file="sigint_popMeans_sumT2drought.txt", header=T)

invDn2dList <- subset(int2dsummary, InvUp==FALSE)
invDn2dList.1 <- invDn2dList$Contig

#identify genes of interest
# #list of sig/not siq
int_invDn2dList <- factor(as.integer(CdifNames %in% invDn2dList.1))
names(int_invDn2dList) <- CdifNames
# str(int_invDn2dList)
# summary(int_invDn2dList)

MFinvDn2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=int_invDn2dList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
MFinvDn2dGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 49  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 48  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635 
MFinvDn2d_resultpc <- runTest(MFinvDn2dGOdata, algorithm = "parentchild", statistic = "fisher")
MFinvDn2d_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 6 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 48 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 130 
MFinvDn2d_pcRes <- GenTable(MFinvDn2dGOdata, parentchild = MFinvDn2d_resultpc, topNodes=6)
MFinvDn2d_pcRes
write.table(MFinvDn2d_pcRes, file="GOresults_MFsigint_invDn2drought_pc.txt", sep="\t")

# invDn2d_pcRes.more<- GenTable(invDn2dGOdata, parentchild = resultpc, topNodes=50)

# get more info on specific terms, here, top 3 
mget(MFinvDn2d_pcRes[,1], GOTERM)

####MF int sig, up in inv, tmpt 2 control####
int2csummary <- read.table(file="sigint_popMeans_sumT2Control.txt", header=T)

invUp2cList <- subset(int2csummary, InvUp==TRUE)
invUp2cList.1 <- invUp2cList$Contig

#identify genes of interest
# #list of sig/not siq
int_invUp2cList <- factor(as.integer(CdifNames %in% invUp2cList.1))
names(int_invUp2cList) <- CdifNames
# str(int_invUp2cList)
# summary(int_invUp2cList)

MFinvUp2cGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=int_invUp2cList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
MFinvUp2cGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 60  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 59  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635 
MFinvUp2c_resultpc <- runTest(MFinvUp2cGOdata, algorithm = "parentchild", statistic = "fisher")
MFinvUp2c_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 5 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 59 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 128 
MFinvUp2c_pcRes <- GenTable(MFinvUp2cGOdata, parentchild = MFinvUp2c_resultpc, topNodes=5)
MFinvUp2c_pcRes
write.table(MFinvUp2c_pcRes, file="GOresults_MFsigint_invUp2control_pc.txt", sep="\t")

# invUp2c_pcRes.more <- GenTable(invUp2cGOdata, parentchild = resultpc, topNodes=30)
mget(MFinvUp2c_pcRes[,1], GOTERM)
# 
####MF int sig, down in inv, tmpt2 control####
# int2csummary <- read.table(file="sigint_popMeans_sumT2Control.txt", header=T)

invDn2cList <- subset(int2csummary, InvUp==FALSE)
invDn2cList.1 <- invDn2cList$Contig

#identify genes of interest
# #list of sig/not siq
int_invDn2cList <- factor(as.integer(CdifNames %in% invDn2cList.1))
names(int_invDn2cList) <- CdifNames
# str(int_invDn2cList)
# summary(int_invDn2cList)

MFinvDn2cGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=int_invDn2cList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
MFinvDn2cGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 135  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 128  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635 
MFinvDn2c_resultpc <- runTest(MFinvDn2cGOdata, algorithm = "parentchild", statistic = "fisher")
MFinvDn2c_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 10 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 128 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 183
MFinvDn2c_pcRes <- GenTable(MFinvDn2cGOdata, parentchild = MFinvDn2c_resultpc, topNodes=10)
MFinvDn2c_pcRes
write.table(MFinvDn2c_pcRes, file="GOresults_MFsigint_invDn2control_pc.txt", sep="\t")

# invDn2c_pcRes.more <- GenTable(invDn2cGOdata, parentchild = resultpc, topNodes=75)
mget(MFinvDn2c_pcRes[,1], GOTERM) 

####CC int sig####
# intqList <- PC1q$intQ
# names(intqList) <- PC1q$Contig
# head(intqList)
#make topGOdata object, MF
CCintqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                    ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                    allGenes=intqList, #factor describing which genes are of interest/sig, which are not
                    geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                    annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                    nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                    gene2GO=GOmap) #our gene->GO term mapping file

CCintqGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 61024 available genes (all genes from the array):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  0.7364 0.761989 0.5792765 0.78036662 0.471659746  ...
# - 227  significant genes. 
# 
# 33904 feasible genes (genes that can be used in the analysis):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  0.7364 0.761989 0.5792765 0.78036662 0.471659746  ...
# - 178  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 418 
# - number of edges = 937 

#enrichment
CCintresultpc <- runTest(CCintqGOdata, algorithm = "parentchild", statistic = "fisher")
CCintresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 418 GO terms scored: 34 terms with p < 0.01
# Annotation data:
#   Annotated genes: 33904 
# Significant genes: 178 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 178 

(CCintpcRes <- GenTable(CCintqGOdata, parentchild = CCintresultpc, topNodes=34))
# CCintpcRes
write.table(CCintpcRes, file="GOresults_CCsigint_pc.txt", sep="\t")

(CCintpcRes.more <- GenTable(CCintqGOdata, parentchild = CCintresultpc, topNodes=40))

####cc int sig, up in inv, tmpt 2, drought####
# int2dsummary <- read.table(file="sigint_popMeans_sumT2drought.txt", header=T)
# summary(int2dsummary)
# invUp2dList <- subset(int2dsummary, InvUp==TRUE)
# invUp2dList.1 <- invUp2dList$Contig
# 
# #identify genes of interest
# # #list of sig/not siq
# int_invUp2dList <- factor(as.integer(CdifNames %in% invUp2dList.1))
# names(int_invUp2dList) <- CdifNames
# str(int_invUp2dList)
# summary(int_invUp2dList)

CCinvUp2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=int_invUp2dList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
CCinvUp2dGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 146  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 141  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006 
CCinvUp2d_resultpc <- runTest(CCinvUp2dGOdata, algorithm = "parentchild", statistic = "fisher")
CCinvUp2d_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 35 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 141 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 175 

CCinvUp2d_pcRes <- GenTable(CCinvUp2dGOdata, parentchild = CCinvUp2d_resultpc, topNodes=35)
CCinvUp2d_pcRes
write.table(CCinvUp2d_pcRes, file="GOresults_CCsigint_invUp2drought_pc.txt", sep="\t")

# CCinvUp2d_pcRes.more <- GenTable(CCinvUp2dGOdata, parentchild = CCinvUp2d_resultpc, topNodes=70)

# get more info on specific terms, here, top 3 
mget(CCinvUp2d_pcRes[,1], GOTERM)



####cc int sig, down in inv, tmpt 2, drought####
# int2dsummary <- read.table(file="sigint_popMeans_sumT2drought.txt", header=T)
# 
# invDn2dList <- subset(int2dsummary, InvUp==FALSE)
# invDn2dList.1 <- invDn2dList$Contig
# 
# #identify genes of interest
# # #list of sig/not siq
# int_invDn2dList <- factor(as.integer(CdifNames %in% invDn2dList.1))
# names(int_invDn2dList) <- CdifNames
# str(int_invDn2dList)
# summary(int_invDn2dList)

CCinvDn2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=int_invDn2dList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
CCinvDn2dGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 49  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 37  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006
CCinvDn2d_resultpc <- runTest(CCinvDn2dGOdata, algorithm = "parentchild", statistic = "fisher")
CCinvDn2d_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 0 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 37 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 85 

#no sig GO terms!!!
CCinvDn2d_pcRes <- GenTable(CCinvDn2dGOdata, parentchild = CCinvDn2d_resultpc, topNodes=10)
CCinvDn2d_pcRes
# write.table(CCinvDn2d_pcRes, file="GOresults_CCsigint_invDn2drought_pc.txt", sep="\t")

# invDn2d_pcRes.more<- GenTable(invDn2dGOdata, parentchild = resultpc, topNodes=50)

# get more info on specific terms, here, top 3 
# mget(CCinvDn2d_pcRes[1:3,1], GOTERM)


####cc int sig, up in inv, tmpt 2 control####
# int2csummary <- read.table(file="sigint_popMeans_sumT2Control.txt", header=T)
# 
# invUp2cList <- subset(int2csummary, InvUp==TRUE)
# invUp2cList.1 <- invUp2cList$Contig
# 
# #identify genes of interest
# # #list of sig/not siq
# int_invUp2cList <- factor(as.integer(CdifNames %in% invUp2cList.1))
# names(int_invUp2cList) <- CdifNames
# str(int_invUp2cList)
# summary(int_invUp2cList)

CCinvUp2cGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=int_invUp2cList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
CCinvUp2cGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 60  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 48  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006 

CCinvUp2c_resultpc <- runTest(CCinvUp2cGOdata, algorithm = "parentchild", statistic = "fisher")
CCinvUp2c_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 3 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 48 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 110 
CCinvUp2c_pcRes <- GenTable(CCinvUp2cGOdata, parentchild = CCinvUp2c_resultpc, topNodes=3)
CCinvUp2c_pcRes
write.table(CCinvUp2c_pcRes, file="GOresults_CCsigint_invUp2control_pc.txt", sep="\t")

# invUp2c_pcRes.more <- GenTable(invUp2cGOdata, parentchild = resultpc, topNodes=30)

# get more info on specific terms, here, top 3 
mget(CCinvUp2c_pcRes[,1], GOTERM)


####cc int sig, down in inv, tmpt2 control####
# int2csummary <- read.table(file="sigint_popMeans_sumT2Control.txt", header=T)
# 
# invDn2cList <- subset(int2csummary, InvUp==FALSE)
# invDn2cList.1 <- invDn2cList$Contig
# 
# #identify genes of interest
# # #list of sig/not siq
# int_invDn2cList <- factor(as.integer(CdifNames %in% invDn2cList.1))
# names(int_invDn2cList) <- CdifNames
# str(int_invDn2cList)
# summary(int_invDn2cList)
# 
CCinvDn2cGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=int_invDn2cList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
CCinvDn2cGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 135  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 130  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006 

CCinvDn2c_resultpc <- runTest(CCinvDn2cGOdata, algorithm = "parentchild", statistic = "fisher")
CCinvDn2c_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 33 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 130 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 169 
CCinvDn2c_pcRes <- GenTable(CCinvDn2cGOdata, parentchild = CCinvDn2c_resultpc, topNodes=33)
CCinvDn2c_pcRes
write.table(CCinvDn2c_pcRes, file="GOresults_CCsigint_invDn2control_pc.txt", sep="\t")
# 
# invDn2c_pcRes.more <- GenTable(invDn2cGOdata, parentchild = resultpc, topNodes=75)
# get more info on specific terms, here, top 3 
mget(CCinvDn2c_pcRes[,1], GOTERM)

####scanning CC go terms####
# invDn2d_pcRes.more
subset(CCinvDn2d_pcRes, GO.ID %in% CCintpcRes$GO.ID)
# subset(invDn2d_pcRes.more, GO.ID %in% intpcRes$GO.ID)
# invUp2d_pcRes.more
subset(CCinvUp2d_pcRes, GO.ID %in% CCintpcRes$GO.ID)
# subset(invUp2d_pcRes.more, GO.ID %in% intpcRes$GO.ID)

# invDn2c_pcRes.more
subset(CCinvDn2c_pcRes, GO.ID %in% CCintpcRes$GO.ID)
# subset(intpcRes, GO.ID %in% invDn2c_pcRes$GO.ID)
# subset(invDn2c_pcRes.more, GO.ID %in% intpcRes$GO.ID)
# invUp2c_pcRes.more
subset(CCinvUp2c_pcRes, GO.ID %in% CCintpcRes$GO.ID)
# subset(invUp2c_pcRes.more, GO.ID %in% intpcRes$GO.ID)

#####Nearly equal! sig int, tmpt 2 dr, inv>nat?################
int2dsummaryNE <- read.table(file="sigint_popMeans_sumT2droughtNE.txt", header=T)
#for InvDefDn, FALSE means down

####Nearly equal! int sig, up in inv rel to nat, tmpt 2, drought####
NEOrUp2dList <- subset(int2dsummaryNE, InvDefUp==TRUE)
NEOrUp2dList.1 <- NEOrUp2dList$Contig

#identify genes of interest
# #list of sig/not siq
NEOrint_invUp2dList <- factor(as.integer(CdifNames %in% NEOrUp2dList.1))
names(NEOrint_invUp2dList) <- CdifNames

#BP
NEOrUp2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=NEOrint_invUp2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
NEOrUp2dGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 124  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 116  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 
#   
NEOrresultpc <- runTest(NEOrUp2dGOdata, algorithm = "parentchild", statistic = "fisher")
NEOrresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 23 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 116 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 466 

NEOrUp2d_pcRes <- GenTable(NEOrUp2dGOdata, parentchild = NEOrresultpc, topNodes=23)
write.table(NEOrUp2d_pcRes, file="GOresults_sigint_invUpreltonat_T2droughtNE_pc.txt", sep="\t")

(NEOrUp2d_pcRes.more <- GenTable(NEOrUp2dGOdata, parentchild = NEOrresultpc, topNodes=60))#44 if alpha=0.05
subset(NEOrUp2d_pcRes.more, GO.ID %in% c("GO:0006970"))

#MF
NEOrUp2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=NEOrint_invUp2dList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
NEOrUp2dGOdata_mf
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 124  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 117  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635  

NEOrUp2d_resultpc_mf <- runTest(NEOrUp2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
NEOrUp2d_resultpc_mf
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 8 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 117 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 184 

NEOrUp2d_mf_pcRes <- GenTable(NEOrUp2dGOdata_mf, parentchild = NEOrUp2d_resultpc_mf, topNodes=8)
NEOrUp2d_mf_pcRes
write.table(NEOrUp2d_mf_pcRes, file="GOresults_sigint_invUpreltonat_T2droughtMFNE_pc.txt", sep="\t")

(NEOrUp2d_mf_pcRes.more <- GenTable(NEOrUp2dGOdata_mf, parentchild = NEOrUp2d_resultpc_mf, topNodes=20))#14 if alpha=0.05

#cc
NEOrUp2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=NEOrint_invUp2dList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
NEOrUp2dGOdata_cc
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 124  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 119  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006

NEOrUp2d_cc_resultpc <- runTest(NEOrUp2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
NEOrUp2d_cc_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 29 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 119 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 163 

NEOrUp2d_cc_pcRes <- GenTable(NEOrUp2dGOdata_cc, parentchild = NEOrUp2d_cc_resultpc, topNodes=29)
NEOrUp2d_cc_pcRes
write.table(NEOrUp2d_cc_pcRes, file="GOresults_sigint_invUpreltonat_T2droughtCCNE_pc.txt", sep="\t")

(NEOrUp2d_cc_pcRes.more <- GenTable(NEOrUp2dGOdata_cc, parentchild = NEOrUp2d_cc_resultpc, topNodes=50))#43 if alpha=0.05

subset(NEOrUp2d_cc_pcRes.more, GO.ID %in% c("GO:0031090","GO:0005740","GO:0044428","GO:0043229","GO:0044444"))


####Nearly Equal! int sig, down in inv rel to nat, tmpt 2, drought####
NEOrDn2dList <- subset(int2dsummaryNE, InvDefDn==FALSE)#for InvDefDn, FALSE means down
NEOrDn2dList.1 <- NEOrDn2dList$Contig

#identify genes of interest
# #list of sig/not siq
NEOrint_invDn2dList <- factor(as.integer(CdifNames %in% NEOrDn2dList.1))
names(NEOrint_invDn2dList) <- CdifNames
 
NEOrDn2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=NEOrint_invDn2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
NEOrDn2dGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 38  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 34  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862    
  
NEOrDnresultpc <- runTest(NEOrDn2dGOdata, algorithm = "parentchild", statistic = "fisher")
NEOrDnresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 8 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 34 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 263 

(NEOrDn2d_pcRes <- GenTable(NEOrDn2dGOdata, parentchild = NEOrDnresultpc, topNodes=8))
write.table(NEOrDn2d_pcRes, file="GOresults_sigint_invDnreltonat_T2droughtNE_pc.txt", sep="\t")

(NEOrDn2d_pcRes.more<- GenTable(NEOrDn2dGOdata, parentchild = NEOrDnresultpc, topNodes=70))#23 if alpha=0.05
subset(NEOrDn2d_pcRes.more, GO.ID %in% c("GO:0006970"))

#MF
NEOrDn2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=NEOrint_invDn2dList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
NEOrDn2dGOdata_mf
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 38  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 37  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635

NEOrDn2d_mf_resultpc <- runTest(NEOrDn2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
NEOrDn2d_mf_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 5 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 37 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 108 

NEOrDn2d_mf_pcRes <- GenTable(NEOrDn2dGOdata_mf, parentchild = NEOrDn2d_mf_resultpc, topNodes=5)
NEOrDn2d_mf_pcRes
write.table(NEOrDn2d_pcRes, file="GOresults_sigint_invDnreltonat_T2droughtNEmf_pc.txt", sep="\t")

(NEOrDn2d_mf_pcRes.more<- GenTable(NEOrDn2dGOdata_mf, parentchild = NEOrDn2d_mf_resultpc, topNodes=25))#12 if alpha=0.05

#cc
NEOrDn2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                       ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                       allGenes=NEOrint_invDn2dList, #factor describing which genes are of interest/sig, which are not
                       #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                       annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                       nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                       gene2GO=GOmap) #our gene->GO term mapping file
NEOrDn2dGOdata_cc
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 38  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 28  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006

NEOrDn2d_cc_resultpc <- runTest(NEOrDn2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
NEOrDn2d_cc_resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 0 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 28 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 75 
#0 sig terms!!!!
NEOrDn2d_cc_pcRes <- GenTable(NEOrDn2dGOdata_cc, parentchild = NEOrDn2d_cc_resultpc, topNodes=4)
NEOrDn2d_cc_pcRes
write.table(NEOrDn2d_cc_pcRes, file="GOresults_sigint_invDnreltonat_T2droughtNEcc_pc.txt", sep="\t")

(NEOrDn2d_cc_pcRes.more<- GenTable(NEOrDn2dGOdata_cc, parentchild = NEOrDn2d_cc_resultpc, topNodes=25))#4 if alpha=0.05

####compare nearly equal sig int results####
#GO terms enriched for rapidly evolving genes in invasive diffusa
# invsyn <- read.table("GOterm_InvSyn.txt", header=T, sep="\t")

#BP - nearly equal
NEOrUp2d_pcRes$GO.ID %in% pcRes$GO.ID
NEOrUp2d_pcRes.more$GO.ID %in% pcRes$GO.ID
subset(NEOrUp2d_pcRes.more, GO.ID %in% pcRes$GO.ID)

NEOrDn2d_pcRes$GO.ID %in% pcRes$GO.ID
NEOrDn2d_pcRes.more$GO.ID %in% pcRes$GO.ID
subset(NEOrDn2d_pcRes.more, GO.ID %in% pcRes$GO.ID)

#CC - nearly equal
NEOrUp2d_cc_pcRes$GO.ID %in% CCintpcRes$GO.ID
NEOrUp2d_cc_pcRes.more$GO.ID %in% CCintpcRes$GO.ID
subset(NEOrUp2d_cc_pcRes.more, GO.ID %in% CCintpcRes$GO.ID)

NEOrDn2d_cc_pcRes$GO.ID %in% CCintpcRes$GO.ID
NEOrDn2d_cc_pcRes.more$GO.ID %in% CCintpcRes$GO.ID
subset(NEOrDn2d_cc_pcRes.more, GO.ID %in% CCintpcRes$GO.ID)

#MF - nearly equal
NEOrUp2d_mf_pcRes$GO.ID %in% MFpcRes$GO.ID
NEOrUp2d_mf_pcRes.more$GO.ID %in% MFpcRes$GO.ID
subset(NEOrUp2d_mf_pcRes.more, GO.ID %in% MFpcRes$GO.ID)

NEOrDn2d_mf_pcRes$GO.ID %in% MFpcRes$GO.ID
NEOrDn2d_mf_pcRes.more$GO.ID %in% MFpcRes$GO.ID
subset(NEOrDn2d_mf_pcRes.more, GO.ID %in% MFpcRes$GO.ID)



##############################################
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
------------------------- topGOdata object -------------------------
  
  Description:
  -  GO analysis of Cdif microarrays; genes with sig Origin effect 

Ontology:
  -  BP 

61024 available genes (all genes from the array):
  - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
- score :  1 0.191291 1 0.62836337 1  ...
- 587  significant genes. 

37373 feasible genes (genes that can be used in the analysis):
  - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
- score :  1 0.191291 1 0.62836337 1  ...
- 369  significant genes. 

GO graph (nodes with at least  10  genes):
  - a graph with directed edges
- number of nodes = 2483 
- number of edges = 5282 

------------------------- topGOdata object -------------------------

#enrichment
#kay used fisher exact test...
# resultclas <- runTest(oqGOdata, algorithm = "classic", statistic = "fisher")
# resultelim <- runTest(oqGOdata, algorithm = "elim", statistic = "fisher")
# resultwt <- runTest(oqGOdata, algorithm = "weight", statistic = "fisher")
# resultwt01 <- runTest(oqGOdata, algorithm = "weight01", statistic = "fisher")
# resultlea <- runTest(oqGOdata, algorithm = "lea", statistic = "fisher")
resultpc <- runTest(oqGOdata, algorithm = "parentchild", statistic = "fisher")

allResO <- GenTable(oqGOdata, classic = resultclas, elim=resultelim, weight=resultwt, 
                   weight01=resultwt01, lea=resultlea, parentchild=resultpc,
                   orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj
write.table(allResO, file="GOresults_sigOrigin.txt", sep="\t")

# Used parentchild in InvSyn paper, so...
resultpc
Description: GO analysis of Cdif microarrays; genes with sig Origin effect 
Ontology: BP 
'parentchild' algorithm with the 'fisher : joinFun = union' test
2483 GO terms scored: 15 terms with p < 0.01
Annotation data:
  Annotated genes: 37373 
Significant genes: 369 
Min. no. of genes annotated to a GO: 10 
Nontrivial nodes: 1088 


pcRes <- GenTable(oqGOdata, parentchild = resultpc, topNodes=15)
write.table(pcRes, file="GOresults_sigOrigin_pc.txt", sep="\t")

# get more info on specific terms, here, top 3 
mget(pcRes[1:3,1], GOTERM)

####O sig, up in inv####
sigOT0summary <- read.table(file="sigOrigin_popMeans_sumT0.txt", header=T)

invUpoT0List <- subset(sigOT0summary, InvUp==TRUE)
invUpoT0List.1 <- invUpoT0List$Contig

#identify genes of interest
# #list of sig/not siq
int_invUpoT0List <- factor(as.integer(CdifNames %in% invUpoT0List.1))
names(int_invUpoT0List) <- CdifNames
# str(int_invUpoT0List)
# summary(int_invUpoT0List)

invUpoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_invUpoT0List, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
invUpoT0GOdata

#   
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 181  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 173  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 


resultpc <- runTest(invUpoT0GOdata, algorithm = "parentchild", statistic = "fisher")
resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 29 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 173 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 747 

invUpoT0_pcRes_bp <- GenTable(invUpoT0GOdata, parentchild = resultpc, topNodes=29)
write.table(invUpoT0_pcRes_bp, file="GOresults_sigorigin_invUpT0_pc.txt", sep="\t")

invUpoT0_pcRes.more <- GenTable(invUpoT0GOdata, parentchild = resultpc, topNodes=100)

# get more info on specific terms, here, top 3 
# mget(pcRes[1:3,1], GOTERM)

####o sig, down in inv####
sigOT0summary <- read.table(file="sigOrigin_popMeans_sumT0.txt", header=T)

invDnoT0List <- subset(sigOT0summary, InvUp==FALSE)
invDnoT0List.1 <- invDnoT0List$Contig

#identify genes of interest
# #list of sig/not siq
int_invDnoT0List <- factor(as.integer(CdifNames %in% invDnoT0List.1))
names(int_invDnoT0List) <- CdifNames
str(int_invDnoT0List)
summary(int_invDnoT0List)

invDnoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                      ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                      allGenes=int_invDnoT0List, #factor describing which genes are of interest/sig, which are not
                      #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                      annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                      nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                      gene2GO=GOmap) #our gene->GO term mapping file
invDnoT0GOdata
------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 212  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 196  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 

------------------------- topGOdata object -------------------------
  
resultpc <- runTest(invDnoT0GOdata, algorithm = "parentchild", statistic = "fisher")
resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 13 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 196 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 830 

invDnoT0_pcRes <- GenTable(invDnoT0GOdata, parentchild = resultpc, topNodes=13)
write.table(invDnoT0_pcRes, file="GOresults_sigorigin_invDnT0_pc.txt", sep="\t")

invDnoT0_pcRes.more <- GenTable(invDnoT0GOdata, parentchild = resultpc, topNodes=70)

# get more info on specific terms, here, top 3 
# mget(pcRes[1:3,1], GOTERM)


####MF Origin sig####
# oqList <- PC1q$originQ
# names(oqList) <- PC1q$Contig

#make topGOdata object
MFoqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin effect",
                ontology="MF", #MF, CC?
                allGenes=oqList, #factor describing which genes are of interest/sig, which are not
                geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                gene2GO=GOmap) #our gene->GO term mapping file

MFoqGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin effect 
# 
# Ontology:
#   -  MF 
# 
# 61024 available genes (all genes from the array):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  1 0.191291 1 0.62836337 1  ...
# - 587  significant genes. 
# 
# 37190 feasible genes (genes that can be used in the analysis):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  1 0.191291 1 0.62836337 1  ...
# - 365  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1100 
# - number of edges = 1376 

#enrichment
MForesultpc <- runTest(MFoqGOdata, algorithm = "parentchild", statistic = "fisher")
MForesultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1100 GO terms scored: 5 terms with p < 0.01
# Annotation data:
#   Annotated genes: 37190 
# Significant genes: 365 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 390 

MFopcRes <- GenTable(MFoqGOdata, parentchild = MForesultpc, topNodes=5)
MFopcRes
write.table(MFopcRes, file="GOresults_MFsigOrigin_pc.txt", sep="\t")

mget(MFopcRes[,1], GOTERM)


####MF O sig, up in inv####
sigOT0summary <- read.table(file="sigOrigin_popMeans_sumT0.txt", header=T)
invUpoT0List <- subset(sigOT0summary, InvUp==TRUE)
invUpoT0List.1 <- invUpoT0List$Contig

#identify genes of interest
# #list of sig/not siq
int_invUpoT0List <- factor(as.integer(CdifNames %in% invUpoT0List.1))
names(int_invUpoT0List) <- CdifNames
# str(int_invUpoT0List)
# summary(int_invUpoT0List)

MFinvUpoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                      ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                      allGenes=int_invUpoT0List, #factor describing which genes are of interest/sig, which are not
                      #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                      annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                      nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                      gene2GO=GOmap) #our gene->GO term mapping file
MFinvUpoT0GOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 181  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 170  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635 
  
MFinvUporesultpc <- runTest(MFinvUpoT0GOdata, algorithm = "parentchild", statistic = "fisher")
MFinvUporesultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 8 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 170 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 260 

MFinvUpopcRes <- GenTable(MFinvUpoT0GOdata, parentchild = MFinvUporesultpc, topNodes=8)
MFinvUpopcRes
write.table(MFinvUpopcRes, file="GOresults_MFsigorigin_invUpT0_pc.txt", sep="\t")

# invUpoT0_pcRes.more <- GenTable(invUpoT0GOdata, parentchild = resultpc, topNodes=100)

mget(MFinvUpopcRes[,1], GOTERM)

####MF o sig, down in inv####
# sigOT0summary <- read.table(file="sigOrigin_popMeans_sumT0.txt", header=T)
invDnoT0List <- subset(sigOT0summary, InvUp==FALSE)
invDnoT0List.1 <- invDnoT0List$Contig

#identify genes of interest
# #list of sig/not siq
int_invDnoT0List <- factor(as.integer(CdifNames %in% invDnoT0List.1))
names(int_invDnoT0List) <- CdifNames
# str(int_invDnoT0List)
# summary(int_invDnoT0List)

MFinvDnoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                      ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                      allGenes=int_invDnoT0List, #factor describing which genes are of interest/sig, which are not
                      #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                      annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                      nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                      gene2GO=GOmap) #our gene->GO term mapping file
MFinvDnoT0GOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 212  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 193  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635 
MFinvDnoresultpc <- runTest(MFinvDnoT0GOdata, algorithm = "parentchild", statistic = "fisher")
MFinvDnoresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 3 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 193 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 282  

MFinvDno_pcRes <- GenTable(MFinvDnoT0GOdata, parentchild = MFinvDnoresultpc, topNodes=3)
MFinvDno_pcRes
write.table(MFinvDno_pcRes, file="GOresults_MFsigorigin_invDnT0_pc.txt", sep="\t")

# invDnoT0_pcRes.more <- GenTable(invDnoT0GOdata, parentchild = resultpc, topNodes=70)

# get more info on specific terms, here, top 3 
mget(MFinvDno_pcRes[,1], GOTERM)



####CC Origin sig####
# oqList <- PC1q$originQ
# names(oqList) <- PC1q$Contig

#make topGOdata object
CCoqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin effect",
                  ontology="CC", #MF, CC?
                  allGenes=oqList, #factor describing which genes are of interest/sig, which are not
                  geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                  annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                  nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                  gene2GO=GOmap) #our gene->GO term mapping file

CCoqGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin effect 
# 
# Ontology:
#   -  CC 
# 
# 61024 available genes (all genes from the array):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  1 0.191291 1 0.62836337 1  ...
# - 587  significant genes. 
# 
# 33904 feasible genes (genes that can be used in the analysis):
#   - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
# - score :  1 0.191291 1 0.62836337 1  ...
# - 340  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 418 
# - number of edges = 937 

#enrichment
CCoresultpc <- runTest(CCoqGOdata, algorithm = "parentchild", statistic = "fisher")
CCoresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 418 GO terms scored: 5 terms with p < 0.01
# Annotation data:
#   Annotated genes: 33904 
# Significant genes: 340 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 226

CCopcRes <- GenTable(CCoqGOdata, parentchild = CCoresultpc, topNodes=5)
CCopcRes
write.table(CCopcRes, file="GOresults_CCsigOrigin_pc.txt", sep="\t")

mget(CCopcRes[,1], GOTERM)

####CC O sig, up in inv####
# sigOT0summary <- read.table(file="sigOrigin_popMeans_sumT0.txt", header=T)
# invUpoT0List <- subset(sigOT0summary, InvUp==TRUE)
# invUpoT0List.1 <- invUpoT0List$Contig

#identify genes of interest
# #list of sig/not siq
# int_invUpoT0List <- factor(as.integer(CdifNames %in% invUpoT0List.1))
# names(int_invUpoT0List) <- CdifNames
# str(int_invUpoT0List)
# summary(int_invUpoT0List)

CCinvUpoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=int_invUpoT0List, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
CCinvUpoT0GOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 181  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 164  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006 

CCinvUporesultpc <- runTest(CCinvUpoT0GOdata, algorithm = "parentchild", statistic = "fisher")
CCinvUporesultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 5 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 164 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 170

CCinvUpopcRes <- GenTable(CCinvUpoT0GOdata, parentchild = CCinvUporesultpc, topNodes=5)
CCinvUpopcRes
write.table(CCinvUpopcRes, file="GOresults_CCsigorigin_invUpT0_pc.txt", sep="\t")

CCinvUpopcRes.more <- GenTable(CCinvUpoT0GOdata, parentchild = CCinvUporesultpc, topNodes=25)
CCinvUpopcRes.more

mget(CCinvUpopcRes[,1], GOTERM)

####CC o sig, down in inv####
# # sigOT0summary <- read.table(file="sigOrigin_popMeans_sumT0.txt", header=T)
# invDnoT0List <- subset(sigOT0summary, InvUp==FALSE)
# invDnoT0List.1 <- invDnoT0List$Contig
# 
# #identify genes of interest
# # #list of sig/not siq
# int_invDnoT0List <- factor(as.integer(CdifNames %in% invDnoT0List.1))
# names(int_invDnoT0List) <- CdifNames
# # str(int_invDnoT0List)
# # summary(int_invDnoT0List)

CCinvDnoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=int_invDnoT0List, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
CCinvDnoT0GOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 212  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 175  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006 

CCinvDnoresultpc <- runTest(CCinvDnoT0GOdata, algorithm = "parentchild", statistic = "fisher")
CCinvDnoresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 10 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 175 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 191 

CCinvDno_pcRes <- GenTable(CCinvDnoT0GOdata, parentchild = CCinvDnoresultpc, topNodes=10)
CCinvDno_pcRes
write.table(CCinvDno_pcRes, file="GOresults_CCsigorigin_invDnT0_pc.txt", sep="\t")

CCinvDno_pcRes.more <- GenTable(CCinvDnoT0GOdata, parentchild = CCinvDnoresultpc, topNodes=50)
CCinvDno_pcRes.more

# get more info on specific terms, here, top 3 
mget(CCinvDno_pcRes[,1], GOTERM)

####Nearly equal! O sig, up in inv####
sigOT0summaryNE <- read.table(file="sigOrigin_popMeans_sumT0NE.txt", header=T)

invUpoT0List <- subset(sigOT0summaryNE, InvDefUp==TRUE)
invUpoT0List.1 <- invUpoT0List$Contig

#identify genes of interest
# #list of sig/not siq
int_invUpoT0List <- factor(as.integer(CdifNames %in% invUpoT0List.1))
names(int_invUpoT0List) <- CdifNames

invUpoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                      ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                      allGenes=int_invUpoT0List, #factor describing which genes are of interest/sig, which are not
                      #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                      annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                      nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                      gene2GO=GOmap) #our gene->GO term mapping file
invUpoT0GOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 168  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 161  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 

oresultpc <- runTest(invUpoT0GOdata, algorithm = "parentchild", statistic = "fisher")
oresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 30 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 161 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 687 

invUpoT0_pcRes <- GenTable(invUpoT0GOdata, parentchild = oresultpc, topNodes=30)
write.table(invUpoT0_pcRes, file="GOresults_sigorigin_invUpT0NE_pc.txt", sep="\t")

(invUpoT0_pcRes.more <- GenTable(invUpoT0GOdata, parentchild = oresultpc, topNodes=100)) #67 if alpah=0.05

#MF
invUpoT0GOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=int_invUpoT0List, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
invUpoT0GOdata_mf
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 168  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 159  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635

invUporesultpc_mf <- runTest(invUpoT0GOdata_mf, algorithm = "parentchild", statistic = "fisher")
invUporesultpc_mf
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 7 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 159 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 250 

invUpopcRes_mf <- GenTable(invUpoT0GOdata_mf, parentchild = invUporesultpc_mf, topNodes=7)
write.table(invUpopcRes_mf, file="GOresults_MFsigorigin_invUpT0NE_pc.txt", sep="\t")

(invUpopcRes_mf.more <- GenTable(invUpoT0GOdata_mf, parentchild = invUporesultpc_mf, topNodes=100))#21 if alpha = 0.05

#CC
invUpoT0GOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=int_invUpoT0List, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
invUpoT0GOdata_cc
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 168  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 152  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006 

invUporesultpc_cc <- runTest(invUpoT0GOdata_cc, algorithm = "parentchild", statistic = "fisher")
invUporesultpc_cc
# Description: GO
 analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 5 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 152 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 159

invUpopcRes_cc <- GenTable(invUpoT0GOdata_cc, parentchild = invUporesultpc_cc, topNodes=5)
invUpopcRes_cc
write.table(invUpopcRes_cc, file="GOresults_CCsigorigin_invUpT0NE_pc.txt", sep="\t")

(invUpopcRes_cc.more <- GenTable(invUpoT0GOdata_cc, parentchild = invUporesultpc_cc, topNodes=100))#12 if alpha = 0.05

####Nearly equal! O sig, dn in inv####
sigOT0summaryNE <- read.table(file="sigOrigin_popMeans_sumT0NE.txt", header=T)
invDnoT0List <- subset(sigOT0summaryNE, InvDefDn==FALSE)
invDnoT0List.1 <- invDnoT0List$Contig

#identify genes of interest
# #list of sig/not siq
int_invDnoT0List <- factor(as.integer(CdifNames %in% invDnoT0List.1))
names(int_invDnoT0List) <- CdifNames


invDnoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin effect",
                      ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                      allGenes=int_invDnoT0List, #factor describing which genes are of interest/sig, which are not
                      #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                      annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                      nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                      gene2GO=GOmap) #our gene->GO term mapping file
invDnoT0GOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 205  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 190  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 
# 
#   
#   resultpc <- runTest(invDnoT0GOdata, algorithm = "parentchild", statistic = "fisher")
# resultpc
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: BP 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 2748 GO terms scored: 13 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 56957 
# # Significant genes: 196 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 830 
# 
# invDnoT0_pcRes <- GenTable(invDnoT0GOdata, parentchild = resultpc, topNodes=13)
# write.table(invDnoT0_pcRes, file="GOresults_sigorigin_invDnT0_pc.txt", sep="\t")
# 
# invDnoT0_pcRes.more <- GenTable(invDnoT0GOdata, parentchild = resultpc, topNodes=70)
# 
# #MF
# MFinvDnoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                         ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                         allGenes=int_invDnoT0List, #factor describing which genes are of interest/sig, which are not
#                         #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                         annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                         nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                         gene2GO=GOmap) #our gene->GO term mapping file
# MFinvDnoT0GOdata
# # Description:
# #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # 
# # Ontology:
# #   -  MF 
# # 
# # 60863 available genes (all genes from the array):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 212  significant genes. 
# # 
# # 56630 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 193  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 1305 
# # - number of edges = 1635 
# MFinvDnoresultpc <- runTest(MFinvDnoT0GOdata, algorithm = "parentchild", statistic = "fisher")
# MFinvDnoresultpc
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: MF 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 1305 GO terms scored: 3 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 56630 
# # Significant genes: 193 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 282  
# 
# MFinvDno_pcRes <- GenTable(MFinvDnoT0GOdata, parentchild = MFinvDnoresultpc, topNodes=3)
# MFinvDno_pcRes
# write.table(MFinvDno_pcRes, file="GOresults_MFsigorigin_invDnT0_pc.txt", sep="\t")
# 
# # invDnoT0_pcRes.more <- GenTable(invDnoT0GOdata, parentchild = resultpc, topNodes=70)
# 
# #CC
# CCinvDnoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                         ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                         allGenes=int_invDnoT0List, #factor describing which genes are of interest/sig, which are not
#                         #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                         annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                         nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                         gene2GO=GOmap) #our gene->GO term mapping file
# CCinvDnoT0GOdata
# # Description:
# #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # 
# # Ontology:
# #   -  CC 
# # 
# # 60863 available genes (all genes from the array):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 212  significant genes. 
# # 
# # 51626 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 175  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 452 
# # - number of edges = 1006 
# 
# CCinvDnoresultpc <- runTest(CCinvDnoT0GOdata, algorithm = "parentchild", statistic = "fisher")
# CCinvDnoresultpc
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: CC 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 452 GO terms scored: 10 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 51626 
# # Significant genes: 175 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 191 
# 
# CCinvDno_pcRes <- GenTable(CCinvDnoT0GOdata, parentchild = CCinvDnoresultpc, topNodes=10)
# CCinvDno_pcRes
# write.table(CCinvDno_pcRes, file="GOresults_CCsigorigin_invDnT0_pc.txt", sep="\t")
# 
# CCinvDno_pcRes.more <- GenTable(CCinvDnoT0GOdata, parentchild = CCinvDnoresultpc, topNodes=50)

####compare nearly equal sig origin results####
#remake full int sig data (pcRes, MFpcRes, CCintpcRes) from Cdif_expr_topgo.R if needed to retrieve full GO term info
#for tables only, load:
opcRes <- read.table(file="GOresults_sigOrigin_pc.txt", header=T)
opcRes_mf <- read.table(file="GOresults_MFsigOrigin_pc.txt", header=T)
opcRes_cc <- read.table(file="GOresults_CCsigOrigin_pc.txt", header=T)


##################################################
####scanning GO terms####
# OpcRes <- read.table(file="GOresults_sigint_pc.txt", header=T)
#doesn't seem like you can read from table?  go back and remake
intpcRes <- pcRes
opcRes <- pcRes
invUp2d_pcRes
invDn2d_pcRes
invUpoT0_pcRes
invDnoT0_pcRes

invDn2d_pcRes.more
subset(invDn2d_pcRes, GO.ID %in% intpcRes$GO.ID)
subset(invDn2d_pcRes.more, GO.ID %in% intpcRes$GO.ID)
invUp2d_pcRes.more
subset(invUp2d_pcRes, GO.ID %in% intpcRes$GO.ID)
subset(invUp2d_pcRes.more, GO.ID %in% intpcRes$GO.ID)

invDn2c_pcRes.more
subset(invDn2c_pcRes, GO.ID %in% intpcRes$GO.ID)
# subset(intpcRes, GO.ID %in% invDn2c_pcRes$GO.ID)
subset(invDn2c_pcRes.more, GO.ID %in% intpcRes$GO.ID)
invUp2c_pcRes.more
subset(invUp2c_pcRes, GO.ID %in% intpcRes$GO.ID)
subset(invUp2c_pcRes.more, GO.ID %in% intpcRes$GO.ID)

invUpoT0_pcRes.more
subset(invUpoT0_pcRes.more, GO.ID %in% opcRes$GO.ID)
invDnoT0_pcRes.more
subset(invDnoT0_pcRes.more, GO.ID %in% opcRes$GO.ID)

# get more info on specific terms, here, top 3 
mget(intpcRes[,1], GOTERM)

mget(invUp2d_pcRes[,1], GOTERM)
mget(invDn2d_pcRes[,1], GOTERM)

mget(opcRes[,1], GOTERM)

mget(invUpoT0_pcRes[,1], GOTERM)
mget(invDnoT0_pcRes[,1], GOTERM)

####compare to InvSyn results####
#GO terms enriched for rapidly evolving genes in invasive diffusa
invsyn <- read.table("GOterm_InvSyn.txt", header=T, sep="\t")

intpcRes$GO.ID %in% invsyn$GO.ID
opcRes$GO.ID %in% invsyn$GO.ID
invUp2d_pcRes$GO.ID %in% invsyn$GO.ID
invDn2d_pcRes$GO.ID %in% invsyn$GO.ID
invUpoT0_pcRes$GO.ID %in% invsyn$GO.ID
invDnoT0_pcRes$GO.ID %in% invsyn$GO.ID

PC1pcRes <- read.table(file="GOresults_sigPC1_pc.txt", header=T)
PC1pcRes$GO.ID %in% invsyn$GO.ID

MFpcRes$GO.ID %in% invsyn$GO.ID
CCintpcRes$GO.ID %in% invsyn$GO.ID
MFopcRes$GO.ID %in% invsyn$GO.ID
CCopcRes$GO.ID %in% invsyn$GO.ID





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
------------------------- topGOdata object -------------------------
  
  Description:
  -  GO analysis of Cdif microarrays; genes with sig Trt effect 

Ontology:
  -  BP 

61024 available genes (all genes from the array):
  - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
- score :  0.30369 0.0795214 0.00130195 0.46292497 0.481873948  ...
- 9617  significant genes. 

37373 feasible genes (genes that can be used in the analysis):
  - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
- score :  0.30369 0.0795214 0.00130195 0.46292497 0.481873948  ...
- 6480  significant genes. 

GO graph (nodes with at least  10  genes):
  - a graph with directed edges
- number of nodes = 2483 
- number of edges = 5282 

------------------------- topGOdata object -------------------------
  
#enrichment
#kay used fisher exact test...
# resultclas <- runTest(trtqGOdata, algorithm = "classic", statistic = "fisher")
# resultelim <- runTest(trtqGOdata, algorithm = "elim", statistic = "fisher")
# resultwt <- runTest(trtqGOdata, algorithm = "weight", statistic = "fisher")
# resultwt01 <- runTest(trtqGOdata, algorithm = "weight01", statistic = "fisher")
# resultlea <- runTest(trtqGOdata, algorithm = "lea", statistic = "fisher")
resultpc <- runTest(trtqGOdata, algorithm = "parentchild", statistic = "fisher")

allResTrt <- GenTable(trtqGOdata, classic = resultclas, elim=resultelim, weight=resultwt, 
                   weight01=resultwt01, lea=resultlea, parentchild=resultpc,
                   orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj
write.table(allResTrt, file="GOresults_sigTrt.txt", sep="\t")

# Used parentchild in InvSyn paper, so...
resultpc
Description: GO analysis of Cdif microarrays; genes with sig Trt effect 
Ontology: BP 
'parentchild' algorithm with the 'fisher : joinFun = union' test
2483 GO terms scored: 236 terms with p < 0.01
Annotation data:
  Annotated genes: 37373 
Significant genes: 6480 
Min. no. of genes annotated to a GO: 10 
Nontrivial nodes: 2359 


pcRes <- GenTable(trtqGOdata, parentchild = resultpc, topNodes=236)
write.table(pcRes, file="GOresults_sigTrt_pc.txt", sep="\t")

# get more info on specific terms, here, top 3 
mget(pcRes[1:3,1], GOTERM)

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
------------------------- topGOdata object -------------------------
  
  Description:
  -  GO analysis of Cdif microarrays; genes with sig PC1 effect 

Ontology:
  -  BP 

61024 available genes (all genes from the array):
  - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
- score :  1 0.262873 1 0.31465626 1  ...
- 1111  significant genes. 

37373 feasible genes (genes that can be used in the analysis):
  - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
- score :  1 0.262873 1 0.31465626 1  ...
- 734  significant genes. 

GO graph (nodes with at least  10  genes):
  - a graph with directed edges
- number of nodes = 2483 
- number of edges = 5282 

------------------------- topGOdata object -------------------------

#enrichment
#kay used fisher exact test...
# resultclas <- runTest(pctqGOdata, algorithm = "classic", statistic = "fisher")
# resultelim <- runTest(pcqGOdata, algorithm = "elim", statistic = "fisher")
# resultwt <- runTest(pcqGOdata, algorithm = "weight", statistic = "fisher")
# resultwt01 <- runTest(pcqGOdata, algorithm = "weight01", statistic = "fisher")
# resultlea <- runTest(pcqGOdata, algorithm = "lea", statistic = "fisher")
resultpc <- runTest(pcqGOdata, algorithm = "parentchild", statistic = "fisher")

allResPC1 <- GenTable(pcqGOdata, classic = resultclas, elim=resultelim, weight=resultwt, 
                   weight01=resultwt01, lea=resultlea, parentchild=resultpc,
                   orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj
write.table(allResPC1, file="GOresults_sigPC1.txt", sep="\t")

# Used parentchild in InvSyn paper, so...
resultpc
Description: GO analysis of Cdif microarrays; genes with sig PC1 effect 
Ontology: BP 
'parentchild' algorithm with the 'fisher : joinFun = union' test
2483 GO terms scored: 22 terms with p < 0.01
Annotation data:
  Annotated genes: 37373 
Significant genes: 734 
Min. no. of genes annotated to a GO: 10 
Nontrivial nodes: 1413 

pcRes <- GenTable(pcqGOdata, parentchild = resultpc, topNodes=22)
write.table(pcRes, file="GOresults_sigPC1_pc.txt", sep="\t")

# get more info on specific terms, here, top 3 
mget(pcRes[1:3,1], GOTERM)

####################################
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