#C.diffusa expression - topGO
#comparing inv drought tmpt0 to inv dr tmpt2
#2/17/15
#with help from Kay, see email chain "questions about microarray analysis" starting 11/14/14

#load data
PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") #contigs and q values
#pop means
popintdf <- read.table("PC1_sigint_popMeans.txt", header=T)
popOdf <- read.table("PC1_sigOrigin_popMeans.txt", header=T)


####up expression in sig int genes, in drought plants tmpt 2, relative to dr plants tmpt 0####
library(reshape2)
library(plyr)

popintdfDr <- subset(popintdf, Trt=="drought"&Tmpt%in%c(0,2), select=c(1,2,4,10:236))
head(popintdfDr)
testdir <- reshape(popintdfDr,direction="long", varying=list(ExprVal=c(4:230)), times=colnames(popintdfDr[4:230]))
head(testdir)
testdir5 <- ddply(testdir, .(time,Origin,Tmpt), summarise, mean(Contig1007))
head(testdir5)

testdir5inv <- subset(testdir5, Origin=="inv", select=c(1,3,4))
testdir5nat <- subset(testdir5, Origin=="nat", select=c(1,3,4))

testdir6inv <- reshape(testdir5inv, direction="wide", idvar="time", timevar="Tmpt")
head(testdir6inv)
colnames(testdir6inv)[2] <- "InvExprValT0"
colnames(testdir6inv)[3] <- "InvExprValT2"
colnames(testdir6inv)[1] <- "Contig"
testdir6inv$T2Up <- testdir6inv$InvExprValT2 > testdir6inv$InvExprValT0

intdrInvsummary <- testdir6inv
head(intdrInvsummary)
summary(intdrInvsummary)
write.table(intdrInvsummary, file="sigint_popMeans_sumT2droughtInv.txt", sep="\t")

testdir6nat <- reshape(testdir5nat, direction="wide", idvar="time", timevar="Tmpt")
head(testdir6nat)
colnames(testdir6nat)[2] <- "NatExprValT0"
colnames(testdir6nat)[3] <- "NatExprValT2"
colnames(testdir6nat)[1] <- "Contig"
testdir6nat$T2Up <- testdir6nat$NatExprValT2 > testdir6nat$NatExprValT0

intdrNatsummary <- testdir6nat
head(intdrNatsummary)
summary(intdrNatsummary)
write.table(intdrNatsummary, file="sigint_popMeans_sumT2droughtNat.txt", sep="\t")

####allowing for nearly equal things####
#floating points in comparisons... don't do equals, just think of less than
close_enough <- function(x, y, tolerance=sqrt(.Machine$double.eps)) {
  abs(x - y) <= tolerance
}

intdrInvsum <- intdrInvsummary
intdrInvsum$T2eq <- close_enough(intdrInvsum$InvExprValT2, intdrInvsum$InvExprValT0, tolerance=0.1)
intdrInvsum$T2DefUp <- intdrInvsum$T2Up 
intdrInvsum[intdrInvsum$T2eq==TRUE,]$T2DefUp <- FALSE
intdrInvsum$T2DefDn <- intdrInvsum$T2Up#FALSE means down
intdrInvsum[intdrInvsum$T2eq==TRUE,]$T2DefDn <- TRUE

intdrNatsum <- intdrNatsummary
intdrNatsum$T2eq <- close_enough(intdrNatsum$NatExprValT2, intdrNatsum$NatExprValT0, tolerance=0.1)
intdrNatsum$T2DefUp <- intdrNatsum$T2Up 
intdrNatsum[intdrNatsum$T2eq==TRUE,]$T2DefUp <- FALSE
intdrNatsum$T2DefDn <- intdrNatsum$T2Up#FALSE means down
intdrNatsum[intdrNatsum$T2eq==TRUE,]$T2DefDn <- TRUE

####topGO analysis of subsets####
source("http://bioconductor.org/biocLite.R")
biocLite(c("topGO", "ALL", "genefilter", "multtest"))

library(topGO)
library(qvalue)

#make custom annotation obj
GOmap <- readMappings(file = "~/GitHub/Cdiff_expression/GOanalysis/GOmap_Cendif1_Cendif2.txt_flip")
str(head(GOmap))

#all contigs, i.e. gene universe
CdifNames <- names(GOmap)
head(CdifNames)

#set alpha level
topDiffGenes <- function(allScore) {
  return(allScore < 0.05)
}

#remake full int sig data (pcRes, MFpcRes) from Cdif_expr_topgo.R


####int sig, inv, up in tmpt 2 rel to tmpt 0, drought####
intdrInvsummary <- read.table(file="sigint_popMeans_sumT2droughtInv.txt", header=T)
summary(intdrInvsummary)
invUp2dList <- subset(intdrInvsummary, T2Up==TRUE)
invUp2dList.1 <- invUp2dList$Contig

#identify genes of interest
# #list of sig/not siq
int_invUp2dList <- factor(as.integer(CdifNames %in% invUp2dList.1))
names(int_invUp2dList) <- CdifNames
str(int_invUp2dList)
summary(int_invUp2dList)

#BP
invUp2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_invUp2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
invUp2dGOdata
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 65  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 58  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 
------------------------- topGOdata object -------------------------
  
InvdT2resultpc <- runTest(invUp2dGOdata, algorithm = "parentchild", statistic = "fisher")
InvdT2resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 1 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 58 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 467 

invdT2_pcRes <- GenTable(invUp2dGOdata, parentchild = InvdT2resultpc, topNodes=1)
write.table(invdT2_pcRes, file="GOresults_sigint_invUpT2relT0dr_pc.txt", sep="\t")

invdT2_pcRes.more <- GenTable(invUp2dGOdata, parentchild = InvdT2resultpc, topNodes=19)

# get more info on specific terms, here, top 3 
mget(pcRes[1:3,1], GOTERM)

##CC
invUp2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_invUp2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
invUp2dGOdata_cc
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 65  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 51  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006 
------------------------- topGOdata object -------------------------
  
InvdT2resultpc_cc <- runTest(invUp2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
InvdT2resultpc_cc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 1 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 51 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 125

(invdT2_cc_pcRes <- GenTable(invUp2dGOdata_cc, parentchild = InvdT2resultpc_cc, topNodes=1))
write.table(invdT2_cc_pcRes, file="GOresults_sigint_invUpT2relT0drCC_pc.txt", sep="\t")

(invdT2_cc_pcRes.more <- GenTable(invUp2dGOdata_cc, parentchild = InvdT2resultpc_cc, topNodes=13))

##MF
invUp2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_invUp2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
invUp2dGOdata_mf
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 65  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 62  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635 
------------------------- topGOdata object -------------------------
  
InvdT2resultpc_mf <- runTest(invUp2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
InvdT2resultpc_mf
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 7 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 62 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 156 

(invdT2_mf_pcRes <- GenTable(invUp2dGOdata_mf, parentchild = InvdT2resultpc_mf, topNodes=7))
write.table(invdT2_mf_pcRes, file="GOresults_sigint_invUpT2relT0drMF_pc.txt", sep="\t")

(invdT2_mf_pcRes.more <- GenTable(invUp2dGOdata_mf, parentchild = InvdT2resultpc_mf, topNodes=15))


####int sig, nat, up in tmpt 2 rel to tmpt 0, drought####
intdrNatsummary <- read.table(file="sigint_popMeans_sumT2droughtNat.txt", header=T)
summary(intdrNatsummary)
NatUp2dList <- subset(intdrNatsummary, T2Up==TRUE)
NatUp2dList.1 <- NatUp2dList$Contig

#identify genes of interest
# #list of sig/not siq
int_NatUp2dList <- factor(as.integer(CdifNames %in% NatUp2dList.1))
names(int_NatUp2dList) <- CdifNames
str(int_NatUp2dList)
summary(int_NatUp2dList)

#BP
NatUp2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_NatUp2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
NatUp2dGOdata
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 87  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 78  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862
------------------------- topGOdata object -------------------------
  
  NatdT2resultpc <- runTest(NatUp2dGOdata, algorithm = "parentchild", statistic = "fisher")
NatdT2resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 12 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 78 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 484 

(NatdT2_pcRes <- GenTable(NatUp2dGOdata, parentchild = NatdT2resultpc, topNodes=12))
write.table(NatdT2_pcRes, file="GOresults_sigint_NatUpT2relT0dr_pc.txt", sep="\t")

(NatdT2_pcRes.more <- GenTable(NatUp2dGOdata, parentchild = NatdT2resultpc, topNodes=28))

#CC
NatUp2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_NatUp2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
NatUp2dGOdata_cc
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 87  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 72  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006
------------------------- topGOdata object -------------------------
  
  NatdT2resultpc_cc <- runTest(NatUp2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
NatdT2resultpc_cc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 15 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 72 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 134  

(NatdT2_cc_pcRes <- GenTable(NatUp2dGOdata_cc, parentchild = NatdT2resultpc_cc, topNodes=15))
write.table(NatdT2_cc_pcRes, file="GOresults_sigint_NatUpT2relT0drCC_pc.txt", sep="\t")

(NatdT2_cc_pcRes.more <- GenTable(NatUp2dGOdata_cc, parentchild = NatdT2resultpc_cc, topNodes=26))

#MF
NatUp2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_NatUp2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
NatUp2dGOdata_mf
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 87  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 84  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635

------------------------- topGOdata object -------------------------
  
  NatdT2resultpc_mf <- runTest(NatUp2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
NatdT2resultpc_mf
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 9 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 84 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 182  

(NatdT2_mf_pcRes <- GenTable(NatUp2dGOdata_mf, parentchild = NatdT2resultpc_mf, topNodes=9))
write.table(NatdT2_mf_pcRes, file="GOresults_sigint_NatUpT2relT0drMF_pc.txt", sep="\t")

(NatdT2_mf_pcRes.more <- GenTable(NatUp2dGOdata_mf, parentchild = NatdT2resultpc_mf, topNodes=19))

####int sig, inv, down in tmpt 2 rel to tmpt 0, drought####
intdrInvsummary <- read.table(file="sigint_popMeans_sumT2droughtInv.txt", header=T)
summary(intdrInvsummary)
invDn2dList <- subset(intdrInvsummary, T2Up==FALSE)
invDn2dList.1 <- invDn2dList$Contig

#identify genes of interest
# #list of sig/not siq
int_invDn2dList <- factor(as.integer(CdifNames %in% invDn2dList.1))
names(int_invDn2dList) <- CdifNames
str(int_invDn2dList)
summary(int_invDn2dList)

#BP
invDn2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_invDn2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
invDn2dGOdata
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 130  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 122  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 

------------------------- topGOdata object -------------------------
  
  InvdT2Dnresultpc <- runTest(invDn2dGOdata, algorithm = "parentchild", statistic = "fisher")
InvdT2Dnresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 27 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 122 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 418 

(invdT2Dn_pcRes <- GenTable(invDn2dGOdata, parentchild = InvdT2Dnresultpc, topNodes=27))
write.table(invdT2Dn_pcRes, file="GOresults_sigint_invDnT2relT0dr_pc.txt", sep="\t")

(invdT2Dn_pcRes.more <- GenTable(invDn2dGOdata, parentchild = InvdT2Dnresultpc, topNodes=48))

##CC
invDn2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=int_invDn2dList, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
invDn2dGOdata_cc
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 130  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 127  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006  
------------------------- topGOdata object -------------------------
  
  InvdT2Dnresultpc_cc <- runTest(invDn2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
InvdT2Dnresultpc_cc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 33 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 127 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 151 

(invdT2Dn_cc_pcRes <- GenTable(invDn2dGOdata_cc, parentchild = InvdT2Dnresultpc_cc, topNodes=33))
write.table(invdT2_cc_pcRes, file="GOresults_sigint_invUpT2relT0drCC_pc.txt", sep="\t")

(invdT2Dn_cc_pcRes.more <- GenTable(invDn2dGOdata_cc, parentchild = InvdT2Dnresultpc_cc, topNodes=40))

##MF
invDn2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=int_invDn2dList, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
invDn2dGOdata_mf
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 130  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 125  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635 
------------------------- topGOdata object -------------------------
  
  InvdT2Dnresultpc_mf <- runTest(invDn2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
InvdT2Dnresultpc_mf
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 7 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 125 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 174

(invdT2Dn_mf_pcRes <- GenTable(invDn2dGOdata_mf, parentchild = InvdT2Dnresultpc_mf, topNodes=7))
write.table(invdT2Dn_mf_pcRes, file="GOresults_sigint_invDnT2relT0drMF_pc.txt", sep="\t")

(invdT2Dn_mf_pcRes.more <- GenTable(invDn2dGOdata_mf, parentchild = InvdT2Dnresultpc_mf, topNodes=14))


####int sig, nat, up in tmpt 2 rel to tmpt 0, drought####
intdrNatsummary <- read.table(file="sigint_popMeans_sumT2droughtNat.txt", header=T)
summary(intdrNatsummary)
NatDn2dList <- subset(intdrNatsummary, T2Up==FALSE)
NatDn2dList.1 <- NatDn2dList$Contig

#identify genes of interest
# #list of sig/not siq
int_NatDn2dList <- factor(as.integer(CdifNames %in% NatDn2dList.1))
names(int_NatDn2dList) <- CdifNames
str(int_NatDn2dList)
summary(int_NatDn2dList)

#BP
NatDn2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=int_NatDn2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
NatDn2dGOdata
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 108  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 102  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 

------------------------- topGOdata object -------------------------
  
  NatdT2Dnresultpc <- runTest(NatDn2dGOdata, algorithm = "parentchild", statistic = "fisher")
NatdT2Dnresultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 39 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 102 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 376 

(NatdT2Dn_pcRes <- GenTable(NatDn2dGOdata, parentchild = NatdT2Dnresultpc, topNodes=39))
write.table(NatdT2Dn_pcRes, file="GOresults_sigint_NatDnT2relT0dr_pc.txt", sep="\t")

(NatdT2Dn_pcRes.more <- GenTable(NatDn2dGOdata, parentchild = NatdT2Dnresultpc, topNodes=54))

#CC
NatDn2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=int_NatDn2dList, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
NatDn2dGOdata_cc
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 108  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 106  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006
------------------------- topGOdata object -------------------------
  
  NatdT2Dnresultpc_cc <- runTest(NatDn2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
NatdT2Dnresultpc_cc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 31 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 106 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 142  

(NatdT2Dn_cc_pcRes <- GenTable(NatDn2dGOdata_cc, parentchild = NatdT2Dnresultpc_cc, topNodes=31))
write.table(NatdT2Dn_cc_pcRes, file="GOresults_sigint_NatDnT2relT0drCC_pc.txt", sep="\t")

(NatdT2Dn_cc_pcRes.more <- GenTable(NatDn2dGOdata_cc, parentchild = NatdT2Dnresultpc_cc, topNodes=45))

#MF
NatDn2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=int_NatDn2dList, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
NatDn2dGOdata_mf
------------------------- topGOdata object -------------------------
#   Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 108  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 103  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635 


------------------------- topGOdata object -------------------------
  
  NatdT2Dnresultpc_mf <- runTest(NatDn2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
NatdT2Dnresultpc_mf
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 5 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 103 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 155  

(NatdT2Dn_mf_pcRes <- GenTable(NatDn2dGOdata_mf, parentchild = NatdT2Dnresultpc_mf, topNodes=5))
write.table(NatdT2Dn_mf_pcRes, file="GOresults_sigint_NatDnT2relT0drMF_pc.txt", sep="\t")

(NatdT2Dn_mf_pcRes.more <- GenTable(NatDn2dGOdata_mf, parentchild = NatdT2Dnresultpc_mf, topNodes=12))


####Nearly equal! int sig, inv, up in tmpt 2 rel to tmpt 0, drought####
#see allowing for nearly equal things above
summary(intdrInvsum)
NEinvUp2dList <- subset(intdrInvsum, T2DefUp==TRUE)
NEinvUp2dList.1 <- NEinvUp2dList$Contig

#identify genes of interest
# #list of sig/not siq
NEint_invUp2dList <- factor(as.integer(CdifNames %in% NEinvUp2dList.1))
names(NEint_invUp2dList) <- CdifNames
# str(int_invUp2dList)
# summary(int_invUp2dList)

#BP
NEinvUp2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=NEint_invUp2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
NEinvUp2dGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 55  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 48  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862 
NEInvdT2resultpc <- runTest(NEinvUp2dGOdata, algorithm = "parentchild", statistic = "fisher")
NEInvdT2resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 1 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 48 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 388 

(NEinvdT2_pcRes <- GenTable(NEinvUp2dGOdata, parentchild = NEInvdT2resultpc, topNodes=1))
write.table(NEinvdT2_pcRes, file="GOresults_sigint_invUpT2relT0drNE_pc.txt", sep="\t")

(NEinvdT2_pcRes.more <- GenTable(NEinvUp2dGOdata, parentchild = NEInvdT2resultpc, topNodes=16))

##CC
NEinvUp2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=NEint_invUp2dList, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
NEinvUp2dGOdata_cc
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  CC 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 55  significant genes. 
# 
# 51626 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 41  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 452 
# - number of edges = 1006
  
NEInvdT2resultpc_cc
####Nearly equal!! int sig, nat, up in tmpt 2 rel to tmpt 0, drought####
 <- runTest(NEinvUp2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
NEInvdT2resultpc_cc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: CC 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 452 GO terms scored: 1 terms with p < 0.01
# Annotation data:
#   Annotated genes: 51626 
# Significant genes: 41 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 101

(NEinvdT2_cc_pcRes <- GenTable(NEinvUp2dGOdata_cc, parentchild = NEInvdT2resultpc_cc, topNodes=1))
write.table(NEinvdT2_cc_pcRes, file="GOresults_sigint_invUpT2relT0drCCNE_pc.txt", sep="\t")

(NEinvdT2_cc_pcRes.more <- GenTable(NEinvUp2dGOdata_cc, parentchild = NEInvdT2resultpc_cc, topNodes=4))

##MF
NEinvUp2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                        ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                        allGenes=NEint_invUp2dList, #factor describing which genes are of interest/sig, which are not
                        #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                        annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                        nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                        gene2GO=GOmap) #our gene->GO term mapping file
NEinvUp2dGOdata_mf
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  MF 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 55  significant genes. 
# 
# 56630 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 53  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 1305 
# - number of edges = 1635
  
NEInvdT2resultpc_mf <- runTest(NEinvUp2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
NEInvdT2resultpc_mf
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: MF 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 1305 GO terms scored: 3 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56630 
# Significant genes: 53 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 131 

(NEinvdT2_mf_pcRes <- GenTable(NEinvUp2dGOdata_mf, parentchild = NEInvdT2resultpc_mf, topNodes=3))
write.table(NEinvdT2_mf_pcRes, file="GOresults_sigint_invUpT2relT0drMFNE_pc.txt", sep="\t")

(NEinvdT2_mf_pcRes.more <- GenTable(NEinvUp2dGOdata_mf, parentchild = NEInvdT2resultpc_mf, topNodes=15))
####Nearly equal!! int sig, nat, up in tmpt 2 rel to tmpt 0, drought####
#see allowing for nearly equal things above
summary(intdrNatsum)
NENatUp2dList <- subset(intdrNatsum, T2DefUp==TRUE)
NENatUp2dList.1 <- NENatUp2dList$Contig

#identify genes of interest
# #list of sig/not siq
NEint_NatUp2dList <- factor(as.integer(CdifNames %in% NENatUp2dList.1))
names(NEint_NatUp2dList) <- CdifNames


#BP
NENatUp2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
                     ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
                     allGenes=NEint_NatUp2dList, #factor describing which genes are of interest/sig, which are not
                     #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
                     annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
                     nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
                     gene2GO=GOmap) #our gene->GO term mapping file
NENatUp2dGOdata
# Description:
#   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# 
# Ontology:
#   -  BP 
# 
# 60863 available genes (all genes from the array):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# - 84  significant genes. 
# 
# 56957 feasible genes (genes that can be used in the analysis):
#   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# - 75  significant genes. 
# 
# GO graph (nodes with at least  10  genes):
#   - a graph with directed edges
# - number of nodes = 2748 
# - number of edges = 5862

NENatdT2resultpc <- runTest(NENatUp2dGOdata, algorithm = "parentchild", statistic = "fisher")
NENatdT2resultpc
# Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# Ontology: BP 
# 'parentchild' algorithm with the 'fisher : joinFun = union' test
# 2748 GO terms scored: 12 terms with p < 0.01
# Annotation data:
#   Annotated genes: 56957 
# Significant genes: 75 
# Min. no. of genes annotated to a GO: 10 
# Nontrivial nodes: 479 

#   
#   NatdT2resultpc <- runTest(NatUp2dGOdata, algorithm = "parentchild", statistic = "fisher")
# NatdT2resultpc
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: BP 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 2748 GO terms scored: 12 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 56957 
# # Significant genes: 78 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 484 
# 
# (NatdT2_pcRes <- GenTable(NatUp2dGOdata, parentchild = NatdT2resultpc, topNodes=12))
# write.table(NatdT2_pcRes, file="GOresults_sigint_NatUpT2relT0dr_pc.txt", sep="\t")
# 
# (NatdT2_pcRes.more <- GenTable(NatUp2dGOdata, parentchild = NatdT2resultpc, topNodes=28))
# 
# #CC
# NatUp2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                         ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                         allGenes=int_NatUp2dList, #factor describing which genes are of interest/sig, which are not
#                         #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                         annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                         nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                         gene2GO=GOmap) #our gene->GO term mapping file
# NatUp2dGOdata_cc
# ------------------------- topGOdata object -------------------------
#   #   Description:
#   #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
#   # 
#   # Ontology:
#   #   -  CC 
#   # 
#   # 60863 available genes (all genes from the array):
#   #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
#   # - 87  significant genes. 
#   # 
#   # 51626 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 72  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 452 
# # - number of edges = 1006
# ------------------------- topGOdata object -------------------------
#   
#   NatdT2resultpc_cc <- runTest(NatUp2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
# NatdT2resultpc_cc
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: CC 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 452 GO terms scored: 15 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 51626 
# # Significant genes: 72 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 134  
# 
# (NatdT2_cc_pcRes <- GenTable(NatUp2dGOdata_cc, parentchild = NatdT2resultpc_cc, topNodes=15))
# write.table(NatdT2_cc_pcRes, file="GOresults_sigint_NatUpT2relT0drCC_pc.txt", sep="\t")
# 
# (NatdT2_cc_pcRes.more <- GenTable(NatUp2dGOdata_cc, parentchild = NatdT2resultpc_cc, topNodes=26))
# 
# #MF
# NatUp2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                         ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                         allGenes=int_NatUp2dList, #factor describing which genes are of interest/sig, which are not
#                         #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                         annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                         nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                         gene2GO=GOmap) #our gene->GO term mapping file
# NatUp2dGOdata_mf
# ------------------------- topGOdata object -------------------------
#   #   Description:
#   #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
#   # 
#   # Ontology:
#   #   -  MF 
#   # 
#   # 60863 available genes (all genes from the array):
#   #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
#   # - 87  significant genes. 
#   # 
#   # 56630 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 84  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 1305 
# # - number of edges = 1635
# 
# ------------------------- topGOdata object -------------------------
#   
#   NatdT2resultpc_mf <- runTest(NatUp2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
# NatdT2resultpc_mf
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: MF 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 1305 GO terms scored: 9 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 56630 
# # Significant genes: 84 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 182  
# 
# (NatdT2_mf_pcRes <- GenTable(NatUp2dGOdata_mf, parentchild = NatdT2resultpc_mf, topNodes=9))
# write.table(NatdT2_mf_pcRes, file="GOresults_sigint_NatUpT2relT0drMF_pc.txt", sep="\t")
# 
# (NatdT2_mf_pcRes.more <- GenTable(NatUp2dGOdata_mf, parentchild = NatdT2resultpc_mf, topNodes=19))
# 
# ####int sig, inv, down in tmpt 2 rel to tmpt 0, drought####
# intdrInvsummary <- read.table(file="sigint_popMeans_sumT2droughtInv.txt", header=T)
# summary(intdrInvsummary)
# invDn2dList <- subset(intdrInvsummary, T2Up==FALSE)
# invDn2dList.1 <- invDn2dList$Contig
# 
# #identify genes of interest
# # #list of sig/not siq
# int_invDn2dList <- factor(as.integer(CdifNames %in% invDn2dList.1))
# names(int_invDn2dList) <- CdifNames
# str(int_invDn2dList)
# summary(int_invDn2dList)
# 
# #BP
# invDn2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                      ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                      allGenes=int_invDn2dList, #factor describing which genes are of interest/sig, which are not
#                      #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                      annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                      nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                      gene2GO=GOmap) #our gene->GO term mapping file
# invDn2dGOdata
# ------------------------- topGOdata object -------------------------
#   #   Description:
#   #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
#   # 
#   # Ontology:
#   #   -  BP 
#   # 
#   # 60863 available genes (all genes from the array):
#   #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
#   # - 130  significant genes. 
#   # 
#   # 56957 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# # - 122  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 2748 
# # - number of edges = 5862 
# 
# ------------------------- topGOdata object -------------------------
#   
#   InvdT2Dnresultpc <- runTest(invDn2dGOdata, algorithm = "parentchild", statistic = "fisher")
# InvdT2Dnresultpc
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: BP 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 2748 GO terms scored: 27 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 56957 
# # Significant genes: 122 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 418 
# 
# (invdT2Dn_pcRes <- GenTable(invDn2dGOdata, parentchild = InvdT2Dnresultpc, topNodes=27))
# write.table(invdT2Dn_pcRes, file="GOresults_sigint_invDnT2relT0dr_pc.txt", sep="\t")
# 
# (invdT2Dn_pcRes.more <- GenTable(invDn2dGOdata, parentchild = InvdT2Dnresultpc, topNodes=48))
# 
# ##CC
# invDn2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                         ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                         allGenes=int_invDn2dList, #factor describing which genes are of interest/sig, which are not
#                         #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                         annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                         nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                         gene2GO=GOmap) #our gene->GO term mapping file
# invDn2dGOdata_cc
# ------------------------- topGOdata object -------------------------
#   #   Description:
#   #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
#   # 
#   # Ontology:
#   #   -  CC 
#   # 
#   # 60863 available genes (all genes from the array):
#   #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
#   # - 130  significant genes. 
#   # 
#   # 51626 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 127  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 452 
# # - number of edges = 1006  
# ------------------------- topGOdata object -------------------------
#   
#   InvdT2Dnresultpc_cc <- runTest(invDn2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
# InvdT2Dnresultpc_cc
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: CC 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 452 GO terms scored: 33 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 51626 
# # Significant genes: 127 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 151 
# 
# (invdT2Dn_cc_pcRes <- GenTable(invDn2dGOdata_cc, parentchild = InvdT2Dnresultpc_cc, topNodes=33))
# write.table(invdT2_cc_pcRes, file="GOresults_sigint_invUpT2relT0drCC_pc.txt", sep="\t")
# 
# (invdT2Dn_cc_pcRes.more <- GenTable(invDn2dGOdata_cc, parentchild = InvdT2Dnresultpc_cc, topNodes=40))
# 
# ##MF
# invDn2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                         ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                         allGenes=int_invDn2dList, #factor describing which genes are of interest/sig, which are not
#                         #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                         annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                         nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                         gene2GO=GOmap) #our gene->GO term mapping file
# invDn2dGOdata_mf
# ------------------------- topGOdata object -------------------------
#   #   Description:
#   #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
#   # 
#   # Ontology:
#   #   -  MF 
#   # 
#   # 60863 available genes (all genes from the array):
#   #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
#   # - 130  significant genes. 
#   # 
#   # 56630 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 125  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 1305 
# # - number of edges = 1635 
# ------------------------- topGOdata object -------------------------
#   
#   InvdT2Dnresultpc_mf <- runTest(invDn2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
# InvdT2Dnresultpc_mf
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: MF 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 1305 GO terms scored: 7 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 56630 
# # Significant genes: 125 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 174
# 
# (invdT2Dn_mf_pcRes <- GenTable(invDn2dGOdata_mf, parentchild = InvdT2Dnresultpc_mf, topNodes=7))
# write.table(invdT2Dn_mf_pcRes, file="GOresults_sigint_invDnT2relT0drMF_pc.txt", sep="\t")
# 
# (invdT2Dn_mf_pcRes.more <- GenTable(invDn2dGOdata_mf, parentchild = InvdT2Dnresultpc_mf, topNodes=14))
# 
# 
# ####int sig, nat, up in tmpt 2 rel to tmpt 0, drought####
# intdrNatsummary <- read.table(file="sigint_popMeans_sumT2droughtNat.txt", header=T)
# summary(intdrNatsummary)
# NatDn2dList <- subset(intdrNatsummary, T2Up==FALSE)
# NatDn2dList.1 <- NatDn2dList$Contig
# 
# #identify genes of interest
# # #list of sig/not siq
# int_NatDn2dList <- factor(as.integer(CdifNames %in% NatDn2dList.1))
# names(int_NatDn2dList) <- CdifNames
# str(int_NatDn2dList)
# summary(int_NatDn2dList)
# 
# #BP
# NatDn2dGOdata <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                      ontology="BP", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                      allGenes=int_NatDn2dList, #factor describing which genes are of interest/sig, which are not
#                      #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                      annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                      nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                      gene2GO=GOmap) #our gene->GO term mapping file
# NatDn2dGOdata
# ------------------------- topGOdata object -------------------------
#   #   Description:
#   #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
#   # 
#   # Ontology:
#   #   -  BP 
#   # 
#   # 60863 available genes (all genes from the array):
#   #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
#   # - 108  significant genes. 
#   # 
#   # 56957 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_rep_c59243 Contig8320  ...
# # - 102  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 2748 
# # - number of edges = 5862 
# 
# ------------------------- topGOdata object -------------------------
#   
#   NatdT2Dnresultpc <- runTest(NatDn2dGOdata, algorithm = "parentchild", statistic = "fisher")
# NatdT2Dnresultpc
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: BP 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 2748 GO terms scored: 39 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 56957 
# # Significant genes: 102 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 376 
# 
# (NatdT2Dn_pcRes <- GenTable(NatDn2dGOdata, parentchild = NatdT2Dnresultpc, topNodes=39))
# write.table(NatdT2Dn_pcRes, file="GOresults_sigint_NatDnT2relT0dr_pc.txt", sep="\t")
# 
# (NatdT2Dn_pcRes.more <- GenTable(NatDn2dGOdata, parentchild = NatdT2Dnresultpc, topNodes=54))
# 
# #CC
# NatDn2dGOdata_cc <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                         ontology="CC", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                         allGenes=int_NatDn2dList, #factor describing which genes are of interest/sig, which are not
#                         #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                         annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                         nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                         gene2GO=GOmap) #our gene->GO term mapping file
# NatDn2dGOdata_cc
# ------------------------- topGOdata object -------------------------
#   #   Description:
#   #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
#   # 
#   # Ontology:
#   #   -  CC 
#   # 
#   # 60863 available genes (all genes from the array):
#   #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
#   # - 108  significant genes. 
#   # 
#   # 51626 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 106  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 452 
# # - number of edges = 1006
# ------------------------- topGOdata object -------------------------
#   
#   NatdT2Dnresultpc_cc <- runTest(NatDn2dGOdata_cc, algorithm = "parentchild", statistic = "fisher")
# NatdT2Dnresultpc_cc
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: CC 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 452 GO terms scored: 31 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 51626 
# # Significant genes: 106 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 142  
# 
# (NatdT2Dn_cc_pcRes <- GenTable(NatDn2dGOdata_cc, parentchild = NatdT2Dnresultpc_cc, topNodes=31))
# write.table(NatdT2Dn_cc_pcRes, file="GOresults_sigint_NatDnT2relT0drCC_pc.txt", sep="\t")
# 
# (NatdT2Dn_cc_pcRes.more <- GenTable(NatDn2dGOdata_cc, parentchild = NatdT2Dnresultpc_cc, topNodes=45))
# 
# #MF
# NatDn2dGOdata_mf <- new("topGOdata", description = "GO analysis of Cdif microarrays; genes with sig Origin*Trt effect",
#                         ontology="MF", #i.e. biological processes, MF (molecular function), CC (cellular component)?
#                         allGenes=int_NatDn2dList, #factor describing which genes are of interest/sig, which are not
#                         #                   geneSel = topDiffGenes, #if above factor contains p or q values, how to set alpha level
#                         annot=annFUN.gene2GO, #func that maps gene names to GO terms; use gene2GO since we provide gene->GO term mapping file
#                         nodeSize=10, #to prune smaller GO terms, which may be artifacts; use range 5 - 10
#                         gene2GO=GOmap) #our gene->GO term mapping file
# NatDn2dGOdata_mf
# ------------------------- topGOdata object -------------------------
#   #   Description:
#   #   -  GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
#   # 
#   # Ontology:
#   #   -  MF 
#   # 
#   # 60863 available genes (all genes from the array):
#   #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
#   # - 108  significant genes. 
#   # 
#   # 56630 feasible genes (genes that can be used in the analysis):
# #   - symbol:  DKTR001_white_c11792 DKTR001_white_s64574 DKUS022_white_c32880 DKUS022_white_c43404 DKUS022_white_rep_c59243  ...
# # - 103  significant genes. 
# # 
# # GO graph (nodes with at least  10  genes):
# #   - a graph with directed edges
# # - number of nodes = 1305 
# # - number of edges = 1635 
# 
# 
# ------------------------- topGOdata object -------------------------
#   
#   NatdT2Dnresultpc_mf <- runTest(NatDn2dGOdata_mf, algorithm = "parentchild", statistic = "fisher")
# NatdT2Dnresultpc_mf
# # Description: GO analysis of Cdif microarrays; genes with sig Origin*Trt effect 
# # Ontology: MF 
# # 'parentchild' algorithm with the 'fisher : joinFun = union' test
# # 1305 GO terms scored: 5 terms with p < 0.01
# # Annotation data:
# #   Annotated genes: 56630 
# # Significant genes: 103 
# # Min. no. of genes annotated to a GO: 10 
# # Nontrivial nodes: 155  
# 
# (NatdT2Dn_mf_pcRes <- GenTable(NatDn2dGOdata_mf, parentchild = NatdT2Dnresultpc_mf, topNodes=5))
# write.table(NatdT2Dn_mf_pcRes, file="GOresults_sigint_NatDnT2relT0drMF_pc.txt", sep="\t")
# 
# (NatdT2Dn_mf_pcRes.more <- GenTable(NatDn2dGOdata_mf, parentchild = NatdT2Dnresultpc_mf, topNodes=12))

####compare results####
#GO terms enriched for rapidly evolving genes in invasive diffusa
invsyn <- read.table("GOterm_InvSyn.txt", header=T, sep="\t")

#BP
NatdT2_pcRes$GO.ID %in% pcRes$GO.ID
NatdT2_pcRes.more$GO.ID %in% pcRes$GO.ID
subset(NatdT2_pcRes.more, GO.ID %in% pcRes$GO.ID)
invdT2_pcRes$GO.ID %in% pcRes$GO.ID
invdT2_pcRes.more$GO.ID %in% pcRes$GO.ID
subset(invdT2_pcRes.more, GO.ID %in% pcRes$GO.ID)

NatdT2Dn_pcRes$GO.ID %in% pcRes$GO.ID
NatdT2Dn_pcRes.more$GO.ID %in% pcRes$GO.ID
subset(NatdT2Dn_pcRes.more, GO.ID %in% pcRes$GO.ID)
invdT2Dn_pcRes$GO.ID %in% pcRes$GO.ID
invdT2Dn_pcRes.more$GO.ID %in% pcRes$GO.ID
subset(invdT2Dn_pcRes.more, GO.ID %in% pcRes$GO.ID)

#CC
NatdT2_cc_pcRes$GO.ID %in% CCintpcRes$GO.ID
NatdT2_cc_pcRes.more$GO.ID %in% CCintpcRes$GO.ID
subset(NatdT2_cc_pcRes.more, GO.ID %in% CCintpcRes$GO.ID)
invdT2_cc_pcRes$GO.ID %in% CCintpcRes$GO.ID
invdT2_cc_pcRes.more$GO.ID %in% CCintpcRes$GO.ID
subset(invdT2_cc_pcRes.more, GO.ID %in% CCintpcRes$GO.ID)

NatdT2Dn_cc_pcRes$GO.ID %in% CCintpcRes$GO.ID
NatdT2Dn_cc_pcRes.more$GO.ID %in% CCintpcRes$GO.ID
subset(NatdT2Dn_cc_pcRes.more, GO.ID %in% CCintpcRes$GO.ID)
invdT2Dn_cc_pcRes$GO.ID %in% CCintpcRes$GO.ID
invdT2Dn_cc_pcRes.more$GO.ID %in% CCintpcRes$GO.ID
subset(invdT2Dn_cc_pcRes.more, GO.ID %in% CCintpcRes$GO.ID)

#MF
NatdT2_mf_pcRes$GO.ID %in% MFpcRes$GO.ID
subset(NatdT2_mf_pcRes, GO.ID %in% MFpcRes$GO.ID)
NatdT2_mf_pcRes.more$GO.ID %in% MFpcRes$GO.ID
subset(NatdT2_mf_pcRes.more, GO.ID %in% MFpcRes$GO.ID)
invdT2_mf_pcRes$GO.ID %in% MFpcRes$GO.ID
invdT2_mf_pcRes.more$GO.ID %in% MFpcRes$GO.ID
subset(invdT2_mf_pcRes.more, GO.ID %in% MFpcRes$GO.ID)

NatdT2Dn_mf_pcRes$GO.ID %in% MFpcRes$GO.ID
NatdT2Dn_mf_pcRes.more$GO.ID %in% MFpcRes$GO.ID
subset(NatdT2Dn_mf_pcRes.more, GO.ID %in% MFpcRes$GO.ID)
invdT2Dn_mf_pcRes$GO.ID %in% MFpcRes$GO.ID
invdT2Dn_mf_pcRes.more$GO.ID %in% MFpcRes$GO.ID
subset(invdT2Dn_mf_pcRes.more, GO.ID %in% MFpcRes$GO.ID)

#invsyn
MFpcRes.more$GO.ID %in% invsyn$GO.ID
CCintpcRes.more$GO.ID %in% invsyn$GO.ID
pcRes$GO.ID %in% invsyn$GO.ID
