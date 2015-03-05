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
------------------------- topGOdata object -------------------------
  
  Description:
  -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 

Ontology:
  -  BP 

61024 available genes (all genes from the array):
  - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
- score :  0.7364 0.761989 0.5792765 0.78036662 0.471659746  ...
- 227  significant genes. 

37373 feasible genes (genes that can be used in the analysis):
  - symbol:  Contig1 Contig10 Contig100 Contig1000 Contig10000  ...
- score :  0.7364 0.761989 0.5792765 0.78036662 0.471659746  ...
- 180  significant genes. 

GO graph (nodes with at least  10  genes):
  - a graph with directed edges
- number of nodes = 2483 
- number of edges = 5282 

------------------------- topGOdata object -------------------------

#enrichment
#kay used fisher exact test...
# resultclas <- runTest(intqGOdata, algorithm = "classic", statistic = "fisher")
# resultelim <- runTest(intqGOdata, algorithm = "elim", statistic = "fisher")
# resultwt <- runTest(intqGOdata, algorithm = "weight", statistic = "fisher")
# resultwt01 <- runTest(intqGOdata, algorithm = "weight01", statistic = "fisher")
# resultlea <- runTest(intqGOdata, algorithm = "lea", statistic = "fisher")

resultpc <- runTest(intqGOdata, algorithm = "parentchild", statistic = "fisher")

# 
allRes <- GenTable(intqGOdata, classic = resultclas, elim=resultelim, weight=resultwt, 
                   weight01=resultwt01, lea=resultlea, parentchild=resultpc,
                   orderBy = "weight", ranksOf = "classic", topNodes = 50) #can summarize multiple result obj
write.table(allRes, file="GOresults_sigint.txt", sep="\t")

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
>resultpc

Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
Ontology: BP 
'parentchild' algorithm with the 'fisher : joinFun = union' test
2483 GO terms scored: 33 terms with p < 0.01
Annotation data:
  Annotated genes: 37373 
Significant genes: 180 
Min. no. of genes annotated to a GO: 10 
Nontrivial nodes: 639 

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
# ------------------------- topGOdata object -------------------------
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
# ------------------------- topGOdata object -------------------------
  
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

# get more info on specific terms, here, top 3 
# mget(pcRes[1:3,1], GOTERM)


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
str(int_invUpoT0List)
summary(int_invUpoT0List)

invUpoT0GOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_invUpoT0List, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
invUpoT0GOdata
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

------------------------- topGOdata object -------------------------

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

invUpoT0_pcRes <- GenTable(invUpoT0GOdata, parentchild = resultpc, topNodes=29)
write.table(invUpoT0_pcRes, file="GOresults_sigorigin_invUpT0_pc.txt", sep="\t")

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

# get more info on specific terms, here, top 3 
# mget(pcRes[1:3,1], GOTERM)

####scanning GO terms####
# OpcRes <- read.table(file="GOresults_sigint_pc.txt", header=T)
#doesn't seem like you can read from table?  go back and remake
intpcRes <- pcRes
opcRes <- pcRes
invUp2d_pcRes
invDn2d_pcRes
invUpoT0_pcRes
invDnoT0_pcRes

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

####Origin AND trt sig???####
#identify genes of interest
# #list of sig/not siq
# intqList <- factor(as.integer(PC1q$intQsig))
# names(intqList) <- CdifNames
# str(intqList)

# intqList <- PC1q$intQ
PC1q_OT <- subset (PC1q, originQsig==TRUE&trtQsig==TRUE)
OTqList <- PC1q_OT$originQ
# names(intqList) <- PC1q$Contig
# head(intqList)

# #make topGOdata object
# intqGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                   ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                   allGenes=intqList, #factor describing which genes are of interest/sig, which are not
#                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                   annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                   nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                   gene2GO=GOmap) #our gene->GO term mapping file
# 
# intqGOdata


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