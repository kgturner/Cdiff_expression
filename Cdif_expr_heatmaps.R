#C.diffusa expression - fun with heatmaps
#1/19/15

####look at results####
#local
PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") 
# PC1qdr <- read.table("lme4_qval_PC1_dr.txt", header=T, sep="\t") 
# Latq <- read.table("lme4_qval_lat.txt", header=T, sep="\t")
exprs.df <- read.table("C:/Users/Kat/Documents/GitHub/Cdiff_expression/Cdifexprs_lme4dat.txt", header=T, sep="\t") #takes a long time!
test <- exprs.df[,c(1:20)]
#make data longways?
library(reshape2)
# test2 <- reshape(test, idvar="tagged", direction="long", 
#                varying=list(m.date=c(36,55,16), lfl=c(9,13,22), lfw=c(10,14,23), lfc=c(8,12,26)),
#                v.names=c("m.date","lfl", "lfw","lfc"))

PC1q_intsig <- subset(PC1q,intQsig==TRUE ) #contigs with sig Origin*Trt
PC1q_Osig <- subset(PC1q,originQsig==TRUE&intQsig==FALSE )#contigs with sig Origin
PC1q_sig <- subset(PC1q, covQsig==TRUE)
PC1q_trtsig <- subset(PC1q, trtQsig==TRUE&intQsig==FALSE)#contigs with sig trt
# PC1q_trtsig <- subset(PC1q, trtQsig==TRUE)

pc1intV <- as.vector(PC1q_intsig$Contig)
PC1q_intsigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1intV)
PC1q_intsigdf <-cbind(exprs.df[,1:14], PC1q_intsigdf)

pc1oV <- as.vector(PC1q_Osig$Contig)
PC1q_Osigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1oV)
PC1q_Osigdf <-cbind(exprs.df[,1:14], PC1q_Osigdf)

pc1trtV <- as.vector(PC1q_trtsig$Contig)
#
PC1q_trtsigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1trtV)
PC1q_trtsigdf <-cbind(exprs.df[,1:14], PC1q_trtsigdf)

pc1V <- as.vector(PC1q_sig$Contig)
PC1q_sigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1V)
PC1q_sigdf <-cbind(exprs.df[,1:14], PC1q_sigdf)

write.table(PC1q_intsigdf, file="PC1_sigint_df.txt", sep="\t")
write.table(PC1q_Osigdf, file="PC1_sigOrigin_df.txt", sep="\t")
write.table(PC1q_sigdf, file="PC1_sigPC1_df.txt", sep="\t")
write.table(PC1q_trtsigdf, file="PC1_sigtrt_df.txt", sep="\t")

# #lat, not run
# Latq_intsig <- subset(Latq,intQsig==TRUE )
# Latq_Osig <- subset(Latq,originQsig==TRUE&intQsig==FALSE )
# Latq_sig <- subset(Latq, covQsig==TRUE)
#not working
# Latq_intsigdf <- subset(exprs.df, Contig%in%Latq_intsig$Contig)
# Latq_Osigdf <- subset(exprs.df, Contig%in%Latq_Osig$Contig)
# Latq_sigdf <- subset(exprs.df, Contig%in%Latq_sig$Contig)
# write.table(Latq_intsigdf, file="Lat_sigint_df.txt", sep="\t")
# write.table(Latq_Osigdf, file="Lat_sigOrigin_df.txt", sep="\t")
# write.table(Latq_sigdf, file="Lat_sigPC1_df.txt", sep="\t")


summary(PC1q_Osig)
# summary(Latq_Osig)



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

# makeEset<-function(eSet, annt){
#   #Creating an ExpressionSet from eSet, a normalized gene expression matrix
#   # and annt, a data.frame containing annotation
#   metadata <- data.frame(labelDescription = colnames(annt), row.names=colnames(annt))
#   phenoData<-new("AnnotatedDataFrame", data=annt, varMetadata=metadata)
#   if (inherits(eSet, "data.frame")) eSet= as.matrix(eSet)
#   if (inherits(eSet, "ExpressionSet")) eSet=exprs(eSet)
#   data.eSet<-new("ExpressionSet", exprs=eSet, phenoData=phenoData)
#   print(varLabels(data.eSet))
#   return(data.eSet)
# }
# PC1inteset <- makeEset(PC1q_intM[,15:241],PC1q_intM[,1:14])
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

####scaled last time point only####
intM2 <- as.matrix(t(subset(PC1q_intsigdf, Tmpt==2,select=c(15:241))))
intdes2 <- subset(PC1q_intsigdf, Tmpt==2,select=c(1:14))
intdes2$OriginTrt <- as.factor(paste0(intdes2$Origin, "_", intdes2$Trt))
inteset2 <- new("ExpressionSet", phenoData = as(intdes2, "AnnotatedDataFrame"),exprs = as.matrix(intM2))
inteset2<-inteset2[,order(intdes2$OriginTrt)]
intdes2<-intdes2[order(intdes2$OriginTrt),]
exprs(inteset2) <- exprs(inteset2)[,row.names(intdes2)]
all(colnames(exprs(inteset2))==row.names(intdes2))

inteset2_sc<-inteset2
exprs(inteset2_sc)<- t(scale(t( exprs(inteset2_sc) )))


cols<-as.character(as.integer(intdes2$OriginTrt))


# by.cols<-c("blue", "#0000DD",  "#0000BB", "#000099", "#000077","#000055",
#            "#000033", "black", "#333300","#555500", "#777700","#999900", "#BBBB00","#DDDD00", "yellow")
row_distance = dist(exprs(inteset2_sc), method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
# col_distance = dist(t(exprs(inteset3_sc)), method = "euclidean")
# col_cluster = hclust(col_distance, method = "complete")

par(mar=c(12,12,12,12))
heatmap.2(exprs(inteset2_sc), trace="none", ColSideColors = cols, 
          dendrogram="row", na.color="grey50",
          margin=c(12,9),col= "heat.colors", 
          Rowv=as.dendrogram(row_cluster), Colv=NA,
          symm=TRUE, cexRow = 0.6,
          cexCol = 1.2, key = TRUE, keysize=1.5,lwid=c(1,15), lhei=c(1,4),labCol=intdes2$SampleID,
          main = "Sig Origin*Trt by origin")
dev.off()
#in natives, block pattern shift apparent between control and drought
#invasives a bit of a jumble

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
# col_distance = dist(t(exprs(inteset3_sc)), method = "euclidean")
# col_cluster = hclust(col_distance, method = "complete")

png("Cdifexprs_heatmap_popMeanT2.png", width=800, height=800, pointsize = 12)
par(mar=c(12,12,12,12))
heatmap.2(exprs(popinteset2_sc), trace="none", ColSideColors = cols, 
          dendrogram="row", na.color="grey50",
          margin=c(12,9),col= "heat.colors", 
          Rowv=as.dendrogram(row_cluster), Colv=NA,
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