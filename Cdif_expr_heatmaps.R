#C.diffusa expression - fun with heatmaps
#1/19/15

####heatmaps and clustering####
PC1q_intsigdf <- read.table("PC1_sigint_df.txt", header=T)
PC1q_Osigdf <- read.table("PC1_sigOrigin_df.txt", header=T)
PC1q_sigdf <- read.table("PC1_sigPC1_df.txt", header=T)
PC1q_trtsigdf <- read.table("PC1_sigtrt_df.txt", header=T)
# 
# PC1q_intM <- as.matrix(t(PC1q_intsigdf[,c(15:241)]))
# # levels(grdat$Origin)[levels(grdat$Origin)=="inv"] <- "Invasive"
# # levels(grdat$Origin)[levels(grdat$Origin)=="nat"] <- "Native"
# orgs <- as.character(as.vector(as.numeric(PC1q_intsigdf$Origin)))
# 
# heatmap.2(PC1q_intM,trace="none", ColSideColors=orgs, scale="none")

#using biobase
source("http://bioconductor.org/biocLite.R")
biocLite(c("Biobase"))
library("Biobase")
library("RColorBrewer")
library("gplots")


####scaled last time point, control, sig int only####
intM2 <- as.matrix(t(subset(PC1q_intsigdf, Trt=="control"&Tmpt==2,select=c(15:241))))
intdes2 <- subset(PC1q_intsigdf, Trt=="control"&Tmpt==2,select=c(1:14))
# intdes2$OriginTrt <- as.factor(paste0(intdes2$Origin, "_", intdes2$Trt))
inteset2 <- new("ExpressionSet", phenoData = as(intdes2, "AnnotatedDataFrame"),exprs = as.matrix(intM2))
inteset2<-inteset2[,order(intdes2$Origin, intdes2$Pop)]
intdes2<-intdes2[order(intdes2$Origin, intdes2$Pop),]
exprs(inteset2) <- exprs(inteset2)[,row.names(intdes2)]
all(colnames(exprs(inteset2))==row.names(intdes2))

inteset2_sc<-inteset2
exprs(inteset2_sc)<- t(scale(t( exprs(inteset2_sc) )))


cols<-as.character(as.integer(intdes2$Origin))

my_palette <- colorRampPalette(c("royalblue", "black", "gold"))(n = 299)
col_breaks = c(seq(-2.5,-1,length=100), # for red
               seq(-1,1,length=100), # for yellow
               seq(1,2.5,length=100)) # for green

# by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
#            "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(inteset2_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
col_distance = dist(t(exprs(inteset2_sc)), method = "euclidean")
col_cluster = hclust(col_distance, method = "complete")

pdf("SigInt_control_tmpt2_heatmap.pdf", useDingbats=FALSE)
# par(mar=c(12,12,12,12))
heatmap.2(exprs(inteset2_sc), trace="none", ColSideColors = cols, 
          dendrogram="both", na.color="grey50",
          margin=c(12,9),
          col= my_palette, 
          breaks=col_breaks,
          Rowv=as.dendrogram(row_cluster), Colv=as.dendrogram(col_cluster),
#           symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1,labRow=NA,
#           lwid=c(1,15), lhei=c(1,4),
          labCol=intdes2$PopTrtPool,
          main = "Genes with significant Origin*Trt \nControl, Tmpt 2")
dev.off()


####scaled last time point, drought, sig int only####
intM2 <- as.matrix(t(subset(PC1q_intsigdf, Trt=="drought"&Tmpt==2,select=c(15:241))))
intdes2 <- subset(PC1q_intsigdf, Trt=="drought"&Tmpt==2,select=c(1:14))
# intdes2$OriginTrt <- as.factor(paste0(intdes2$Origin, "_", intdes2$Trt))
inteset2 <- new("ExpressionSet", phenoData = as(intdes2, "AnnotatedDataFrame"),exprs = as.matrix(intM2))
inteset2<-inteset2[,order(intdes2$Origin, intdes2$Pop)]
intdes2<-intdes2[order(intdes2$Origin, intdes2$Pop),]
exprs(inteset2) <- exprs(inteset2)[,row.names(intdes2)]
all(colnames(exprs(inteset2))==row.names(intdes2))

inteset2_sc<-inteset2
exprs(inteset2_sc)<- t(scale(t( exprs(inteset2_sc) )))


cols<-as.character(as.integer(intdes2$Origin))

my_palette <- colorRampPalette(c("royalblue", "black", "gold"))(n = 299)
col_breaks = c(seq(-2.5,-1,length=100), # for red
               seq(-1,1,length=100), # for yellow
               seq(1,2.5,length=100)) # for green


# by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
#            "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(inteset2_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
col_distance = dist(t(exprs(inteset2_sc)), method = "euclidean")
col_cluster = hclust(col_distance, method = "complete")

pdf("SigInt_drought_tmpt2_heatmap.pdf", useDingbats=FALSE)
# par(mar=c(12,12,12,12))
heatmap.2(exprs(inteset2_sc), trace="none", ColSideColors = cols, 
          dendrogram="both", na.color="grey50",
          margin=c(12,9),
          col= my_palette,
          breaks=col_breaks,
          Rowv=as.dendrogram(row_cluster), Colv=as.dendrogram(col_cluster),
#           symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1,labRow=NA,
#           lwid=c(1,15), lhei=c(1,4),
          labCol=intdes2$PopTrtPool,
          main = "Genes with significant Origin*Trt \nDrought, Tmpt 2")
dev.off()

####scaled first time point, both trt equivalent, sig O only####
oM2 <- as.matrix(t(subset(PC1q_Osigdf, Tmpt==0,select=c(15:599))))
odes2 <- subset(PC1q_Osigdf, Tmpt==0,select=c(1:14))
oeset2 <- new("ExpressionSet", phenoData = as(odes2, "AnnotatedDataFrame"),exprs = as.matrix(oM2))
oeset2<-oeset2[,order(odes2$Origin, odes2$Pop)]
odes2<-odes2[order(odes2$Origin, odes2$Pop),]
exprs(oeset2) <- exprs(oeset2)[,row.names(odes2)]
all(colnames(exprs(oeset2))==row.names(odes2))

oeset2_sc<-oeset2
exprs(oeset2_sc)<- t(scale(t( exprs(oeset2_sc) )))

cols<-as.character(as.integer(odes2$Origin))

my_palette <- colorRampPalette(c("royalblue", "black", "gold"))(n = 299)
col_breaks = c(seq(-3,-1,length=100), # for red
               seq(-1,1,length=100), # for yellow
               seq(1,3,length=100)) # for green

row_distance = dist(exprs(oeset2_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
col_distance = dist(t(exprs(oeset2_sc)), method = "euclidean")
col_cluster = hclust(col_distance, method = "complete")

pdf("SigOrigin_bothtrt_tmpt0_heatmap.pdf", useDingbats=FALSE)
# par(mar=c(12,12,12,12))

heatmap.2(exprs(oeset2_sc), trace="none", ColSideColors = cols, 
          dendrogram="both", na.color="grey50",
          margin=c(12,9),
          col= my_palette, 
          breaks=col_breaks, 
          Rowv=as.dendrogram(row_cluster), Colv=as.dendrogram(col_cluster),
#           cexRow = 0.6,
          labRow=NA,
          cexCol = 1.2, 
          key = TRUE, keysize=1,
#           lwid=c(1,15), lhei=c(1,4),
          labCol=odes2$PopTrtPool,
          main = "Genes with significant effect of Origin \nBoth Trt, Tmpt 0")
dev.off()

####scaled last time point, both trt, sig int only####
intM2 <- as.matrix(t(subset(PC1q_intsigdf, Tmpt==2,select=c(15:241))))
intdes2 <- subset(PC1q_intsigdf, Tmpt==2,select=c(1:14))
# intdes2$OriginTrt <- as.factor(paste0(intdes2$Origin, "_", intdes2$Trt))
inteset2 <- new("ExpressionSet", phenoData = as(intdes2, "AnnotatedDataFrame"),exprs = as.matrix(intM2))
inteset2<-inteset2[,order(intdes2$Origin, intdes2$Pop)]
intdes2<-intdes2[order(intdes2$Origin, intdes2$Pop),]
exprs(inteset2) <- exprs(inteset2)[,row.names(intdes2)]
all(colnames(exprs(inteset2))==row.names(intdes2))

inteset2_sc<-inteset2
exprs(inteset2_sc)<- t(scale(t( exprs(inteset2_sc) )))


cols<-as.character(as.integer(intdes2$Origin))

my_palette <- colorRampPalette(c("royalblue", "black", "gold"))(n = 299)
col_breaks = c(seq(-2,-1,length=100), # for red
               seq(-1,1,length=100), # for yellow
               seq(1,2,length=100)) # for green


# # by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
# #            "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(inteset2_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
col_distance = dist(t(exprs(inteset2_sc)), method = "euclidean")
col_cluster = hclust(col_distance, method = "complete")
# 
pdf("SigInt_bothtrt_tmpt2_heatmap.pdf", useDingbats=FALSE)
# # par(mar=c(12,12,12,12))
# par(cex.main=1)
heatmap.2(exprs(inteset2_sc), trace="none", 
          ColSideColors = cols, 
          na.color="grey50", dendrogram="both", 
          margin=c(12,9),
          col= my_palette,
          breaks=col_breaks,
          Rowv=as.dendrogram(row_cluster), Colv=as.dendrogram(col_cluster),
          #           hclustfun=function(x) hclust(x,method="complete"),
          #           distfun=function(x) dist(x,method="euclidean"),
          #           Rowv=TRUE,Colv=TRUE,
#           symm=TRUE, cexRow = 0.6,
          cexCol = 1, key = TRUE, keysize=1,
#           lwid=c(1,15), lhei=c(1,4),
          labCol=intdes2$SampleID,labRow=NA,
          main = "Genes with significant effect of Origin*Trt \nBoth Trt, Tmpt 2")
dev.off()


###############maybe?################
# ####scaled last time point only, split treatments, sig int and sig O genes####
intO <- merge(PC1q_intsigdf, PC1q_Osigdf) #combine sig int and sig Origin genes
intOC <- subset(intO, Trt=="control"&Tmpt==2)
intODr <- subset(intO, Trt=="drought"&Tmpt==2)

#control
intOCM <- as.matrix(t(subset(intOC, select=c(15:826))))
intOCdes <- subset(intOC, select=c(1:14))
# intOCdes$OriginTrt <- as.factor(paste0(intdes2$Origin, "_", intdes2$Trt))
intOCeset <- new("ExpressionSet", phenoData = as(intOCdes, "AnnotatedDataFrame"),exprs = as.matrix(intOCM))
intOCeset<-intOCeset[,order(intOCdes$Origin, intOCdes$Pop)]
intOCdes<-intOCdes[order(intOCdes$Origin,intOCdes$Pop),]
exprs(intOCeset) <- exprs(intOCeset)[,row.names(intOCdes)]
all(colnames(exprs(intOCeset))==row.names(intOCdes))
# 
intOCeset_sc<-intOCeset
exprs(intOCeset_sc)<- t(scale(t( exprs(intOCeset_sc) )))
# 
# 
cols<-as.character(as.integer(intOCdes$Origin))
# 
# 
# # by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
# #            "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(intOCeset_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete") #should this be euclidean too???
col_distance = dist(t(exprs(intOCeset_sc)), method = "euclidean")
col_cluster = hclust(col_distance, method = "complete")
#
pdf("SigOandInt_control_tmpt2_heatmap.pdf", useDingbats=FALSE)
# par(mar=c(2,2,2,2))
heatmap.2(exprs(intOCeset_sc), trace="none", ColSideColors = cols, 
          dendrogram="row", na.color="grey50",
          margin=c(12,9),col= "heat.colors", 
          Rowv=as.dendrogram(row_cluster), Colv=as.dendrogram(col_cluster),
          symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),labCol=intOCdes$SampleID,
          main = "Control, Tmpt 2, Sig Origin*Trt and Sig Origin, sorted by origin")
dev.off()

####scaled last time point, control, sig origin only####
intM2 <- as.matrix(t(subset(PC1q_Osigdf, Trt=="control"&Tmpt==2,select=c(15:241))))
intdes2 <- subset(PC1q_Osigdf, Trt=="control"&Tmpt==2,select=c(1:14))
# intdes2$OriginTrt <- as.factor(paste0(intdes2$Origin, "_", intdes2$Trt))
inteset2 <- new("ExpressionSet", phenoData = as(intdes2, "AnnotatedDataFrame"),exprs = as.matrix(intM2))
inteset2<-inteset2[,order(intdes2$Origin, intdes2$Pop)]
intdes2<-intdes2[order(intdes2$Origin, intdes2$Pop),]
exprs(inteset2) <- exprs(inteset2)[,row.names(intdes2)]
all(colnames(exprs(inteset2))==row.names(intdes2))

inteset2_sc<-inteset2
exprs(inteset2_sc)<- t(scale(t( exprs(inteset2_sc) )))


cols<-as.character(as.integer(intdes2$Origin))


# by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
#            "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(inteset2_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
# col_distance = dist(t(exprs(inteset3_sc)), method = "euclidean")
# col_cluster = hclust(col_distance, method = "complete")
pdf("SigO_control_tmpt2_heatmap.pdf", useDingbats=FALSE)
# par(mar=c(12,12,12,12))
heatmap.2(exprs(inteset2_sc), trace="none", ColSideColors = cols, 
          dendrogram="row", na.color="grey50",
          margin=c(12,9),col= "heat.colors", 
          Rowv=as.dendrogram(row_cluster), Colv=NA,
          symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),labCol=intdes2$PopTrtPool,
          main = "Control, Tmpt 2,Sig Origin by origin")
dev.off()





#######trials below#####
####pc1 int####
PC1q_intM <- as.matrix(t(PC1q_intsigdf[,c(15:241)]))

all(colnames(t(PC1q_intsigdf[,c(15:241)]))==row.names(PC1q_intsigdf))


PC1inteset <- new("ExpressionSet", phenoData = as(PC1q_intsigdf[,1:14], "AnnotatedDataFrame"),exprs = as.matrix(PC1q_intM))
PC1inteset<-PC1inteset[,order( PC1q_intsigdf$Trt, PC1q_intsigdf$Origin )]
PC1q_intsigdf<-PC1q_intsigdf[order(PC1q_intsigdf$Trt),]
exprs(PC1inteset) <- exprs(PC1inteset)[,row.names(PC1q_intsigdf)]
all(colnames(exprs(PC1inteset))==row.names(PC1q_intsigdf))
cols<-as.character(as.integer(PC1q_intsigdf$Trt))

par(mar=c(5,6,4,2)+0.1,mgp=c(3,1,0), oma=c(0,0,0,0))
by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
           "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
heatmap.2(exprs(PC1inteset), trace="none", ColSideColors = cols, 
          dendrogram="none", na.color="grey50",
          margin=c(4,5),col= by.cols, hclustfun = hcf, 
          Rowv=NA, Colv=NA, symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),
          main = "Sig Origin*Trt by trt", scale="none")
dev.off()

# Set up labels and legend
Trt <- levels(PC1q_intsigdf$Trt) 
legend_labelA <- Trt[order(Trt)]
legend_colsA <- as.character(as.integer(order(Trt)))
legend("left",legend=legend_labelA, fill=legend_colsA)

#####control only####
PC1q_intMC <- as.matrix(t(PC1q_intsigdf[PC1q_intsigdf$Trt=="control",c(15:241)]))
PC1q_intdesC <- PC1q_intsigdf[PC1q_intsigdf$Trt=="control",1:14]
PC1intesetC <- new("ExpressionSet", phenoData = as(PC1q_intdesC, "AnnotatedDataFrame"),exprs = as.matrix(PC1q_intMC))
PC1intesetC<-PC1intesetC[,order( PC1q_intdesC$Origin )]
PC1q_intdesC<-PC1q_intdesC[order(PC1q_intdesC$Origin),]
exprs(PC1intesetC) <- exprs(PC1intesetC)[,row.names(PC1q_intdesC)]
all(colnames(exprs(PC1intesetC))==row.names(PC1q_intdesC))
cols<-as.character(as.integer(PC1q_intdesC$Origin))

par(mar=c(5,6,4,2)+0.1,mgp=c(3,1,0), oma=c(0,0,0,0))
by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
           "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
heatmap.2(exprs(PC1intesetC), trace="none", ColSideColors = cols, 
          dendrogram="none", na.color="grey50",
          margin=c(4,5),col= by.cols, hclustfun = hcf, 
          Rowv=NA, Colv=NA, symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),
          main = "Sig Origin*Trt by origin", scale="none")
dev.off()

# # Set up labels and legend
# Trt <- levels(PC1q_intdesC$Trt) 
# legend_labelA <- Trt[order(Trt)]
# legend_colsA <- as.character(as.integer(order(Trt)))
# legend("left",legend=legend_labelA, fill=legend_colsA)

####middle timepoint,  only####

intM1 <- as.matrix(t(subset(PC1q_intsigdf, Tmpt==1,select=c(15:241))))
intdes1 <- subset(PC1q_intsigdf, Tmpt==1,select=c(1:14))
intdes1$OriginTrt <- as.factor(paste0(intdes1$Origin, "_", intdes1$Trt))
inteset1 <- new("ExpressionSet", phenoData = as(intdes1, "AnnotatedDataFrame"),exprs = as.matrix(intM1))
inteset1<-inteset1[,order(intdes1$OriginTrt)]
intdes1<-intdes1[order(intdes1$OriginTrt),]
exprs(inteset1) <- exprs(inteset1)[,row.names(intdes1)]
all(colnames(exprs(inteset1))==row.names(intdes1))
cols<-as.character(as.integer(intdes1$OriginTrt))

par(mar=c(12,12,12,12))
by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
           "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
heatmap.2(exprs(inteset1), trace="none", ColSideColors = cols, 
          dendrogram="none", na.color="grey50",
          margin=c(12,9),col= by.cols, hclustfun = hcf, 
          Rowv=NA, Colv=NA, symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),labCol=intdes1$SampleID,
          main = "Sig Origin*Trt by origin", scale="none")
dev.off()

# # Set up labels and legend
# Trt <- levels(PC1q_intdesC$Trt) 
# legend_labelA <- Trt[order(Trt)]
# legend_colsA <- as.character(as.integer(order(Trt)))
# legend("left",legend=legend_labelA, fill=legend_colsA)

#### scaled middle timepoint,  only####
intM1 <- as.matrix(t(subset(PC1q_intsigdf, Tmpt==1,select=c(15:241))))
intdes1 <- subset(PC1q_intsigdf, Tmpt==1,select=c(1:14))
intdes1$OriginTrt <- as.factor(paste0(intdes1$Origin, "_", intdes1$Trt))
inteset1 <- new("ExpressionSet", phenoData = as(intdes1, "AnnotatedDataFrame"),exprs = as.matrix(intM1))
inteset1<-inteset1[,order(intdes1$OriginTrt)]
intdes1<-intdes1[order(intdes1$OriginTrt),]
exprs(inteset1) <- exprs(inteset1)[,row.names(intdes1)]
all(colnames(exprs(inteset1))==row.names(intdes1))

inteset1_sc<-inteset1
exprs(inteset1_sc)<- t(scale(t( exprs(inteset1_sc) )))


cols<-as.character(as.integer(intdes1$OriginTrt))


by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
           "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(inteset1_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
col_distance = dist(t(exprs(inteset1_sc)), method = "euclidean")
col_cluster = hclust(col_distance, method = "complete")

par(mar=c(12,12,12,12))
heatmap.2(exprs(inteset1_sc), trace="none", ColSideColors = cols, 
          dendrogram="row", na.color="grey50",
          margin=c(12,9),col= "heat.colors", 
          Rowv=as.dendrogram(row_cluster), Colv=NA,
          symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),labCol=intdes1$SampleID,
          main = "Sig Origin*Trt by origin")
dev.off()

# # Set up labels and legend
# Trt <- levels(PC1q_intdesC$Trt) 
# legend_labelA <- Trt[order(Trt)]
# legend_colsA <- as.character(as.integer(order(Trt)))
# legend("left",legend=legend_labelA, fill=legend_colsA)



####scaled first timepoint####
#should  be no difference due to treat
intM0 <- as.matrix(t(subset(PC1q_intsigdf, Tmpt==0,select=c(15:241))))
intdes0 <- subset(PC1q_intsigdf, Tmpt==0,select=c(1:14))
intdes0$OriginTrt <- as.factor(paste0(intdes0$Origin, "_", intdes2$Trt))
inteset0 <- new("ExpressionSet", phenoData = as(intdes0, "AnnotatedDataFrame"),exprs = as.matrix(intM0))
inteset0<-inteset0[,order(intdes0$OriginTrt)]
intdes0<-intdes0[order(intdes0$OriginTrt),]
exprs(inteset0) <- exprs(inteset0)[,row.names(intdes0)]
all(colnames(exprs(inteset0))==row.names(intdes0))

inteset0_sc<-inteset0
exprs(inteset0_sc)<- t(scale(t( exprs(inteset0_sc) )))


cols<-as.character(as.integer(intdes0$OriginTrt))


# by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
#            "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(inteset0_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
# col_distance = dist(t(exprs(inteset3_sc)), method = "euclidean")
# col_cluster = hclust(col_distance, method = "complete")

par(mar=c(12,12,12,12))
heatmap.2(exprs(inteset0_sc), trace="none", ColSideColors = cols, 
          dendrogram="row", na.color="grey50",
          margin=c(12,9),col= "heat.colors", 
          Rowv=as.dendrogram(row_cluster), Colv=NA,
          symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),labCol=intdes0$SampleID,
          main = "Sig Origin*Trt by origin first tmpt")
dev.off()

####average by population####
library(plyr)
test <- PC1q_intsigdf[,c(1:34)]
test2 <- ddply(test, .(Pop, Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))

popintdf <- ddply(PC1q_intsigdf, .(Pop, Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))
popintdf$OriginTrt <- as.factor(paste0(popintdf$Origin, "_", popintdf$Trt))
popintdf$PopTrtTmpt <- as.factor(paste0(popintdf$Pop, "_", popintdf$Trt,"_",popintdf$Tmpt))

write.table(popintdf, file="PC1_sigint_popMeans.txt", sep="\t")

popintdf <- read.table("PC1_sigint_popMeans.txt", header=T)

#for last tmpt
popintM2 <- as.matrix(t(subset(popintdf, Tmpt==2,select=c(10:236))))
popintdes2 <- subset(popintdf, Tmpt==2,select=c(1:9,237:238))
# popintdes2$OriginTrt <- as.factor(paste0(popintdes2$Origin, "_", popintdes2$Trt))
# popintdes2$PopTrtTmpt <- as.factor(paste0(popintdes2$Pop, "_", popintdes2$Trt,"_",popintdes2$Tmpt))

popinteset2 <- new("ExpressionSet", phenoData = as(popintdes2, "AnnotatedDataFrame"),exprs = as.matrix(popintM2))
popinteset2<-popinteset2[,order(popintdes2$OriginTrt)]
popintdes2<-popintdes2[order(popintdes2$OriginTrt),]
exprs(popinteset2) <- exprs(popinteset2)[,row.names(popintdes2)]
all(colnames(exprs(popinteset2))==row.names(popintdes2))

popinteset2_sc<-popinteset2
exprs(popinteset2_sc)<- t(scale(t( exprs(popinteset2_sc) )))


cols<-as.character(as.integer(popintdes2$OriginTrt))


# by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
#            "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(popinteset2_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
col_distance = dist(t(exprs(popinteset2_sc)), method = "euclidean")
col_cluster = hclust(col_distance, method = "complete")

png("Cdifexprs_heatmap_popMeanT2.png", width=800, height=800, pointsize = 12)
par(mar=c(12,12,12,12))
heatmap.2(exprs(popinteset2_sc), trace="none", ColSideColors = cols, 
          dendrogram="row", na.color="grey50",
          margin=c(12,9),col= "heat.colors", 
          Rowv=as.dendrogram(row_cluster), Colv=as.dendrogram(col_cluster),
          symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),labCol=popintdes2$PopTrtTmpt,
          main = "pop means for unigenes w/ Sig Origin*Trt \nLast time point")

dev.off()

#for last time point, inv only
popintM2I <- as.matrix(t(subset(popintdf, Tmpt==2&Origin=="inv",select=c(10:236))))
popintdes2I <- subset(popintdf, Tmpt==2&Origin=="inv",select=c(1:9,237:238))

popinteset2I <- new("ExpressionSet", phenoData = as(popintdes2I, "AnnotatedDataFrame"),exprs = as.matrix(popintM2I))
popinteset2I<-popinteset2I[,order(popintdes2I$OriginTrt)]
popintdes2I<-popintdes2I[order(popintdes2I$OriginTrt),]
exprs(popinteset2I) <- exprs(popinteset2I)[,row.names(popintdes2I)]
all(colnames(exprs(popinteset2I))==row.names(popintdes2I))

popinteset2I_sc<-popinteset2I
exprs(popinteset2I_sc)<- t(scale(t( exprs(popinteset2I_sc) )))


cols<-as.character(as.integer(popintdes2I$OriginTrt))


# by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
#            "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(popinteset2I_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
# col_distance = dist(t(exprs(inteset3_sc)), method = "euclidean")
# col_cluster = hclust(col_distance, method = "complete")

png("Cdifexprs_heatmap_popMeanT2_Inv.png", width=800, height=800, pointsize = 12)
par(mar=c(12,12,12,12))
heatmap.2(exprs(popinteset2I_sc), trace="none", ColSideColors = cols, 
          dendrogram="row", na.color="grey50",
          margin=c(12,9),col= "heat.colors", 
          Rowv=as.dendrogram(row_cluster), Colv=NA,
          symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),labCol=popintdes2I$PopTrtTmpt,
          main = "pop means for unigenes w/ Sig Origin*Trt \nLast time point, Invasives only")

dev.off()


####assorted plots####
library(ggplot2)
library(plyr)
# library(reshape2)

grdat <- PC1q_intsigdf[,c(1:34)]
levels(grdat$Origin)[levels(grdat$Origin)=="inv"] <- "Invasive"
levels(grdat$Origin)[levels(grdat$Origin)=="nat"] <- "Native"
#already longways!
# grdat <- reshape(grdat, idvar="PopTrtPool", direction="long", 
#                varying=list(m.date=c(36,55,16), lfl=c(9,13,22), lfw=c(10,14,23), lfc=c(8,12,26)),
#                v.names=c("m.date","lfl", "lfw","lfc"))
grd_20 <- ddply(grdat, .(Pop, Origin, Trt, PC1,Tmpt), summarize, popCount=length(Pop), 
                popContig1007=mean(Contig1007,na.rm = TRUE),popContig10572=mean(Contig10572,na.rm = TRUE),
                popContig10573=mean(Contig10573,na.rm = TRUE),popContig10722=mean(Contig10722,na.rm = TRUE),
                popContig10774=mean(Contig10774,na.rm = TRUE),popContig10847=mean(Contig10847,na.rm = TRUE),
                popContig11096=mean(Contig11096,na.rm = TRUE),popContig11106=mean(Contig11106,na.rm = TRUE),
                popContig11151=mean(Contig11151,na.rm = TRUE),popContig11376=mean(Contig11376,na.rm = TRUE),
                popContig11541=mean(Contig11541,na.rm = TRUE),popContig11691=mean(Contig11691,na.rm = TRUE),
                popContig11926=mean(Contig11926,na.rm = TRUE),popContig12049=mean(Contig12049,na.rm = TRUE),
                popContig12069=mean(Contig12069,na.rm = TRUE),popContig12082=mean(Contig12082,na.rm = TRUE),
                popContig12462=mean(Contig12462,na.rm = TRUE),popContig12728=mean(Contig12728,na.rm = TRUE),
                popContig12765=mean(Contig12765,na.rm = TRUE),popContig12795=mean(Contig12795,na.rm = TRUE)) #avg per timepoint
#single Interaction
#facet by tmpt?
p1007 <- ggplot(grd_20,aes(Trt, popContig1007,color=Origin))+geom_point(aes(shape=Origin, color=Origin), size=3)+ 
  facet_grid(~ Tmpt,scales="free_y")+
  geom_smooth(method=glm, se=TRUE)+ #ylim(0,1)+
  #coord_cartesian(ylim = c(0, 1.02)) +
  xlab("PC1")+ylab("Population mean popContig1007")+ 
  theme_bw() #+
#   theme(legend.justification=c(0,0), legend.position=c(0,0),
#         legend.title = element_text(size=14, face="bold"),
#         legend.text = element_text(size = 13))
p1007
#avg over time?
grd_20timeless <- ddply(grdat, .(Pop, Origin, Trt, PC1), summarize, popCount=length(Pop), 
                        popContig1007=mean(Contig1007,na.rm = TRUE),popContig10572=mean(Contig10572,na.rm = TRUE),
                        popContig10573=mean(Contig10573,na.rm = TRUE),popContig10722=mean(Contig10722,na.rm = TRUE),
                        popContig10774=mean(Contig10774,na.rm = TRUE),popContig10847=mean(Contig10847,na.rm = TRUE),
                        popContig11096=mean(Contig11096,na.rm = TRUE),popContig11106=mean(Contig11106,na.rm = TRUE),
                        popContig11151=mean(Contig11151,na.rm = TRUE),popContig11376=mean(Contig11376,na.rm = TRUE),
                        popContig11541=mean(Contig11541,na.rm = TRUE),popContig11691=mean(Contig11691,na.rm = TRUE),
                        popContig11926=mean(Contig11926,na.rm = TRUE),popContig12049=mean(Contig12049,na.rm = TRUE),
                        popContig12069=mean(Contig12069,na.rm = TRUE),popContig12082=mean(Contig12082,na.rm = TRUE),
                        popContig12462=mean(Contig12462,na.rm = TRUE),popContig12728=mean(Contig12728,na.rm = TRUE),
                        popContig12765=mean(Contig12765,na.rm = TRUE),popContig12795=mean(Contig12795,na.rm = TRUE)) #avg per timepoint

p1007timeless <- ggplot(grd_20timeless,aes(Trt, popContig1007,color=Origin, group=Origin))+geom_point(aes(shape=Origin, color=Origin), size=3)+ 
  #   facet_grid(~ Tmpt,scales="free_y")+
  geom_smooth(method=glm, se=TRUE)+ #ylim(0,1)+
  #coord_cartesian(ylim = c(0, 1.02)) +
  xlab("PC1")+ylab("Population mean popContig1007")+ 
  theme_bw() #+
#   theme(legend.justification=c(0,0), legend.position=c(0,0),
#         legend.title = element_text(size=14, face="bold"),
#         legend.text = element_text(size = 13))
p1007timeless

#facet by timepoint, PC1 on the x axis, pop means
# #time goes right
# grdat_l <- PC1q_intsig
# levels(grdat_l$Origin)[levels(grdat_l$Origin)=="inv"] <- "Invasive"
# levels(grdat_l$Origin)[levels(grdat_l$Origin)=="nat"] <- "Native"
# grd_l2 <- ddply(grdat_l, .(Pop, Origin, Trt, PC1,time), summarize, popCount=length(Pop), 
#                 poplfc=mean(lfc,na.rm = TRUE)) #avg per timepoint
# pLfc.2 <- ggplot(grd_l2,aes(PC1, poplfc,color=Origin))+geom_point()+ 
#   facet_grid(Trt ~ time,scales="free_y")+
#   geom_smooth(method=glm, se=TRUE)+ #ylim(0,1)+
#   #coord_cartesian(ylim = c(0, 1.02)) +
#   xlab("PC1")+ylab("Population mean leaf number")+ 
#   theme_bw() #+
# #   theme(legend.justification=c(0,0), legend.position=c(0,0),
# #         legend.title = element_text(size=14, face="bold"),
# #         legend.text = element_text(size = 13))
# pLfc.2
#