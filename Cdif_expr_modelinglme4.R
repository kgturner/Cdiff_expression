#Statistical modeling run on cluster (darjeeling) using lme4 package
#11/11/2014

#Install/compile/add to path R version 3.1.1

#set up screen in shell
#screen

#start R session in shell
#R

#install necessary packages
source("http://bioconductor.org/biocLite.R")
#install packages
biocLite(c("limma", "oligo", "pdInfoBuilder", "qvalue"))
#biocLite("limma")
#biocLite("maanova")
#biocLite("pdInfoBuilder")
# biocLite("qvalue")

install.packages("lme4.0", type="both",repos=c("http://lme4.r-forge.r-project.org/repos",getOption("repos")[["CRAN"]]))
#type="both" doesn't work on linux... remove?

#load libraries
library(limma)
library(oligo)
library(pdInfoBuilder)
library(qvalue)
library(lme4.0)


#now install custom package for this experiment. there is no need to load it in after. R will recognize that it's installed once you try to read the xys files
install.packages("~/Centaurea_diffusa_expression/pd.110405.cdiffusa.lz.exp/", repos=NULL, type="source")

####load data####
#read in KNN imputed XYS files
xys.files<-list.xysfiles("~/Centaurea_diffusa_expression/knn_xys", full.names=TRUE)
Cdifxys<-read.xysfiles(xys.files)

#RMA - background correction, normalization, summarization of probes
ppData.Cdifxys=rma(Cdifxys)

#RMA is complete, now extract expression data
exprs.Cdifxys=exprs(ppData.Cdifxys)

# Read hybridization design
design <- read.table("~/Centaurea_diffusa_expression/experimentaldesign.txt", header=T, sep="\t") 
#expt design table includes expression and DNA arrays
"arrayID"  "SampleID" "Pop"      "Origin"   "Trt"      "Exp"      "TrtPool"  "Tmpt"
#arrayID needs to be Array
colnames(design)[1] <- "Array"

exprs.design <- droplevels(subset(design, Exp=="exprs"))

####format data####
#change data matrix to dataframe, and rotate
# test <- as.data.frame(exprs.Cdifxys)
t(head(exprs.Cdifxys))
exprs.df <- as.data.frame(t(exprs.Cdifxys))
row.names(exprs.df) <- gsub("_532\\.xys","",row.names(exprs.df)) 
exprs.df$Array <- as.factor(row.names(exprs.df))

#add climate and lat info to design file
#see Cdif_exprs_climate_pca.R
exprs.clim <- read.table("CdifExprs.BioclimPCAdat.txt", header=TRUE)
exprs.design <- merge(exprs.design, exprs.clim[,c(1,22:27)], all.x=TRUE)
exprs.design$PopTrtPool <- paste0(exprs.design$Pop, "_",exprs.design$TrtPool)
#write design file with climate data
write.table(exprs.design, file="experimentaldesign.txt")

exprs.design <- read.table("C:/Users/Kat/Documents/GitHub/Cdiff_expression/experimentaldesign.txt", header=T, sep="\t") 

#add design info to data frame
exprs.df <- merge(exprs.design,exprs.df, all.y=TRUE)
#write data+design file. quite slow.
write.table(exprs.df, file="Cdifexprs_lme4dat.txt", sep="\t")

# write.table(test, file="test_lme4dat.txt", sep="\t")
# test.1 <- read.table("test_lme4dat.txt", header=T, sep="\t")
####reload data, if neccessary####
#read data table, probably slow
#on cluster
exprs.df <- read.table("~/Centaurea_diffusa_expression/Cdifexprs_lme4dat.txt", header=T, sep="\t") 
#local
exprs.df <- read.table("C:/Users/Kat/Documents/GitHub/Cdiff_expression/Cdifexprs_lme4dat.txt", header=T, sep="\t") 

####loop writing####
# #using lme4.0
# #for loop?
# for (i in colnames(exprs.df[,14:61037])){
#   model1 <- lmer(i ~ Origin*Trt+Latitude+ (Tmpt|SampleID)+(Origin|Pop), family=gaussian, data=exprs.df)
# }

#func?
exprs.LR<- function(trait,df,cov, family=gaussian){
  modeldata<-df[!is.na(df[[trait]]),]
  
  #browser()
  
  model1<-lmer(modeldata[[trait]]  ~ Origin*Trt+modeldata[[cov]]+ (Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)
  model2<-lmer(modeldata[[trait]]  ~ Origin+Trt+ modeldata[[cov]] + (Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)
  model3<-lmer(modeldata[[trait]]  ~ Origin+Trt + (Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)
  model4 <- lmer(modeldata[[trait]]  ~ Trt+ modeldata[[cov]]+(Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=modeldata)

  a1 <- anova(model2,model1) # is interaction sig?
  a2 <- anova(model3,model2) # is covariate sig?
  a3 <- anova(model4, model2) #is origin sig?
  
  pval <- as.data.frame(cbind(Contig=trait,intLRT=a1[[7]][2],covLRT=a2[[7]][2],originLRT=a3[[7]][2]))
  
  return(pval)
}

# test <- exprs.df[, c(1:20)]
# test.LRT <- do.call(rbind,lapply(names(test)[15:20],function(n) exprs.LR(n,df=test,cov="Latitude")))#apply func to all things in list
# test2.LRT <- do.call(rbind,lapply(names(test)[15:20],function(n) exprs.LR(n,df=test)))#apply func to all things in list

# #func testing
# test <- exprs.df[,c(1:20)]
# # test$PopTrtPool <- paste0(test$Pop, "_",test$TrtPool)
# model1<-lmer(Contig1  ~ Origin*Trt+Latitude+(Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=test)
# summary(model1)
# model2<-lmer(Contig1  ~ Origin+Trt+Latitude+(Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=test)
# a1 <- anova(model2,model1)
# model3<-lmer(Contig1  ~ Origin+Trt+ (Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=test)
# a2 <- anova(model3,model2)
# model4 <- lmer(Contig1  ~ Trt+ Latitude+(Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=test)
# a3 <- anova(model4, model2)
# #get out pval of LRT
# a2[[7]][2]
# #paste to new dataframe
# #cols: contig, intLRT, covLRT, originLRT
# #aka: trait, a1[[7]][2], ditto a2, ditto a3
# trait <- "Contig1"
# cov <- "Latitude"
# pval <- as.data.frame(cbind(Contig=trait,intLRT=a1[[7]][2],covLRT=a2[[7]][2],originLRT=a3[[7]][2]))
# 
# model1<-lmer(test[[trait]]  ~ Origin*Trt+test[[cov]]+ (Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=test)

####corrections for multiple tests####
#using:
# test <- exprs.df[, c(1:20)]
# test.LRT <- do.call(rbind,lapply(names(test)[15:20],function(n) exprs.LR(n,df=test,cov="Latitude")))#apply func to all things in list

intP <- as.numeric(as.vector(test.LRT$intLRT))
covP <- as.numeric(as.vector(test.LRT$covLRT))
originP <- as.numeric(as.vector(test.LRT$originLRT))
#from Kay's paper: FDR cuttoff %5 (i.e. fdr.level=0.05)
#otherwise, default from Storey and Tibshirani(2003)
#or maybe use pi0.method="bootstrap" from Storey, Taylor & Siegmund (2004)
intQ <- qvalue(p=intP, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
       smooth.df=3, smooth.log.pi0=FALSE)
covQ <- qvalue(p=covP, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
               smooth.df=3, smooth.log.pi0=FALSE)
test.LRT.Q <- cbind(test.LRT, covQ=covQ$qvalues, covQsig=covQ$significant)

originQ <- qvalue(p=originP, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
               smooth.df=3, smooth.log.pi0=FALSE)


####run models####
#latitude
exprs.LRT.lat <- do.call(rbind,lapply(names(exprs.df)[15:61038],function(n) exprs.LR(n,df=exprs.df,cov="Latitude")))#apply func to all things in list

intQ <- qvalue(p=as.numeric(as.vector(exprs.LRT.lat$intLRT)), lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
               smooth.df=3, smooth.log.pi0=FALSE)
covQ <- qvalue(p=as.numeric(as.vector(exprs.LRT.lat$covLRT)), lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
               smooth.df=3, smooth.log.pi0=FALSE)
originQ <- qvalue(p=as.numeric(as.vector(exprs.LRT.lat$originLRT)), lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
                  smooth.df=3, smooth.log.pi0=FALSE)

exprs.LRT.lat.Q <- cbind(exprs.LRT.lat, intQ=intQ$qvalues, intQsig=intQ$significant,covQ=covQ$qvalues, covQsig=covQ$significant,
                         originQ=originQ$qvalues, originQsig=originQ$significant)
#write pval table. slow? Be sure to include cov name
write.table(exprs.LRT.lat.Q, file="lme4_qval_lat.txt", sep="\t")



#PC1
exprs.LRT.PC1 <- do.call(rbind,lapply(names(exprs.df)[15:61038],function(n) exprs.LR(n,df=exprs.df,cov="PC1")))#apply func to all things in list

intQ <- qvalue(p=as.numeric(as.vector(exprs.LRT.PC1$intLRT)), lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
               smooth.df=3, smooth.log.pi0=FALSE)
covQ <- qvalue(p=as.numeric(as.vector(exprs.LRT.PC1$covLRT)), lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
               smooth.df=3, smooth.log.pi0=FALSE)
originQ <- qvalue(p=as.numeric(as.vector(exprs.LRT.PC1$originLRT)), lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
                  smooth.df=3, smooth.log.pi0=FALSE)

exprs.LRT.PC1.Q <- cbind(exprs.LRT.PC1, intQ=intQ$qvalues, intQsig=intQ$significant,covQ=covQ$qvalues, covQsig=covQ$significant,
                         originQ=originQ$qvalues, originQsig=originQ$significant)
#write pval table. slow? Be sure to include cov name
write.table(exprs.LRT.PC1.Q, file="lme4_qval_PC1.txt", sep="\t")

# #necessary?
# exprs.LRT.PC2 <- do.call(rbind,lapply(names(exprs.df)[15:61038],function(n) exprs.LR(n,df=exprs.df,cov="PC2")))#apply func to all things in list
# #write data+design file. slow? Be sure to include cov name
# write.table(exprs.LRT.pval, file="lme4_pval_PC2.txt")

####look at results####
#local
PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") 
Latq <- read.table("lme4_qval_lat.txt", header=T, sep="\t")
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

pc1intV <- as.vector(PC1q_intsig$Contig)
PC1q_intsigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1intV)
PC1q_intsigdf <-cbind(exprs.df[,1:14], PC1q_intsigdf)
pc1oV <- as.vector(PC1q_Osig$Contig)
PC1q_Osigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1oV)
PC1q_Osigdf <-cbind(exprs.df[,1:14], PC1q_Osigdf)

pc1V <- as.vector(PC1q_sig$Contig)
PC1q_sigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1V)
PC1q_sigdf <-cbind(exprs.df[,1:14], PC1q_sigdf)

write.table(PC1q_intsigdf, file="PC1_sigint_df.txt", sep="\t")
write.table(PC1q_Osigdf, file="PC1_sigOrigin_df.txt", sep="\t")
write.table(PC1q_sigdf, file="PC1_sigPC1_df.txt", sep="\t")

#lat, not run
Latq_intsig <- subset(Latq,intQsig==TRUE )
Latq_Osig <- subset(Latq,originQsig==TRUE&intQsig==FALSE )
Latq_sig <- subset(Latq, covQsig==TRUE)
#not working
# Latq_intsigdf <- subset(exprs.df, Contig%in%Latq_intsig$Contig)
# Latq_Osigdf <- subset(exprs.df, Contig%in%Latq_Osig$Contig)
# Latq_sigdf <- subset(exprs.df, Contig%in%Latq_sig$Contig)
# write.table(Latq_intsigdf, file="Lat_sigint_df.txt", sep="\t")
# write.table(Latq_Osigdf, file="Lat_sigOrigin_df.txt", sep="\t")
# write.table(Latq_sigdf, file="Lat_sigPC1_df.txt", sep="\t")


summary(PC1q_Osig)
summary(Latq_Osig)

####plots####
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

####heatmaps and clustering####
library(gplots)

PC1q_intM <- as.matrix(t(PC1q_intsigdf[,c(15:241)]))
# levels(grdat$Origin)[levels(grdat$Origin)=="inv"] <- "Invasive"
# levels(grdat$Origin)[levels(grdat$Origin)=="nat"] <- "Native"
orgs <- as.character(as.vector(as.numeric(PC1q_intsigdf$Origin)))

heatmap.2(PC1q_intM,trace="none", ColSideColors=orgs, scale="none")

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

#control only
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

# Set up labels and legend
Trt <- levels(PC1q_intdesC$Trt) 
legend_labelA <- Trt[order(Trt)]
legend_colsA <- as.character(as.integer(order(Trt)))
legend("left",legend=legend_labelA, fill=legend_colsA)