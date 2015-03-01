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
#test!
test <- read.table("test_lme4dat.txt", header=T, sep="\t") 

####single contig modeling####
# test <- exprs.df[, c(1:20)]
#test!
test <- read.table("test_lme4dat.txt", header=T, sep="\t") 

modeldata<-test[!is.na(test$Contig100),]
model1<-lmer(Contig100  ~ Origin*Trt+PC1+ (Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=modeldata)
model2<-lmer(Contig100  ~ Origin+Trt+ PC1 + (Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=modeldata)
model3<-lmer(Contig100  ~ Origin+Trt + (Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=modeldata)
model4 <- lmer(Contig100  ~ Trt+ PC1+(Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=modeldata)
modeltrt<-lmer(Contig100  ~ Origin+ PC1 + (Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=modeldata)


(a1 <- anova(model2,model1)) # is interaction sig?
(a2 <- anova(model3,model2)) # is covariate sig?
(a3 <- anova(model4, model2)) #is origin sig?
(a4 <- anova(modeltrt, model2)) #is trt sig?

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
  model4 <- lmer(modeldata[[trait]]  ~ Trt+ modeldata[[cov]]+(Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)

  a1 <- anova(model2,model1) # is interaction sig?
  a2 <- anova(model3,model2) # is covariate sig?
  a3 <- anova(model4, model2) #is origin sig?
  
  pval <- as.data.frame(cbind(Contig=trait,intLRT=a1[[7]][2],covLRT=a2[[7]][2],originLRT=a3[[7]][2]))
  
  return(pval)
}

#test!
test <- read.table("test_lme4dat.txt", header=T, sep="\t") 

test <- exprs.df[, c(1:20)]
test.LRT <- do.call(rbind,lapply(names(test)[15:20],function(n) exprs.LR(n,df=test,cov="Latitude")))#apply func to all things in list
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

exprs.dr.LR<- function(trait,df,cov, family=gaussian){
  modeldata<-df[!is.na(df[[trait]]),]
  
  #browser()
  
  model1<-lmer(modeldata[[trait]]  ~ Origin+Trt+modeldata[[cov]]+ (Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)
  model2<-lmer(modeldata[[trait]]  ~ Origin+ modeldata[[cov]] + (Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)
  
  a1 <- anova(model2,model1) # is Trt sig?
  
  pval <- as.data.frame(cbind(Contig=trait,trtLRT=a1[[7]][2]))
  
  return(pval)
}

#test!
test <- read.table("test_lme4dat.txt", header=T, sep="\t") 
# test <- exprs.df[, c(1:20)]

test.LRT <- do.call(rbind,lapply(names(test)[15:20],function(n) exprs.dr.LR(n,df=test,cov="PC1")))#apply func to all things in list

exprs.pop.LR<- function(trait,df,cov, family=gaussian){
  modeldata<-df[!is.na(df[[trait]]),]
  
  #browser()
  
  model1<-lmer(modeldata[[trait]]  ~ Origin+Trt+modeldata[[cov]]+ (Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)
  model2<-lmer(modeldata[[trait]]  ~ Origin+Trt+ modeldata[[cov]] + (Tmpt|PopTrtPool), family,data=modeldata)
  
  a1 <- anova(model2,model1) # is pop sig?
  
  pval <- as.data.frame(cbind(Contig=trait,popLRT=a1[[7]][2]))
  
  return(pval)
}

#test!
test <- read.table("test_lme4dat.txt", header=T, sep="\t") 

test.LRT <- do.call(rbind,lapply(names(test)[15:20],function(n) exprs.pop.LR(n,df=test,cov="PC1")))#apply func to all things in list

####corrections for multiple tests####
#using:
# test <- exprs.df[, c(1:20)]
test <- read.table("test_lme4dat.txt", header=T, sep="\t") 
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

#trt sig w/ PC1
exprs.dr.LRT.PC1 <- do.call(rbind,lapply(names(exprs.df)[15:61038],function(n) exprs.dr.LR(n,df=exprs.df,cov="PC1")))#apply func to all things in list

trtQ <- qvalue(p=as.numeric(as.vector(exprs.dr.LRT.PC1$trtLRT)), lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
               smooth.df=3, smooth.log.pi0=FALSE)

exprs.dr.LRT.PC1.Q <- cbind(exprs.dr.LRT.PC1, trtQ=trtQ$qvalues, trtQsig=trtQ$significant)
#write pval table. slow? Be sure to include cov name
write.table(exprs.dr.LRT.PC1.Q, file="lme4_qval_PC1_dr.txt", sep="\t")

#add to 1st PC1 table
# test <- cbind(exprs.LRT.PC1.Q, exprs.dr.LRT.PC1.Q)
# write.table(test, file="lme4_qval_PC1.txt", sep="\t")

PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") 
PC1qdr <- read.table("lme4_qval_PC1_dr.txt", header=T, sep="\t") 

test <- merge(PC1q[,1:10], PC1qdr)
write.table(test, file="lme4_qval_PC1.txt", sep="\t")

###
#pop sig w/ PC1
exprs.pop.LRT.PC1 <- do.call(rbind,lapply(names(exprs.df)[15:61038],function(n) exprs.pop.LR(n,df=exprs.df,cov="PC1")))#apply func to all things in list

#start here
popQ <- qvalue(p=as.numeric(as.vector(exprs.pop.LRT.PC1$popLRT)), lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=0.05, robust=FALSE, gui=FALSE, 
               smooth.df=3, smooth.log.pi0=FALSE)

exprs.pop.LRT.PC1.Q <- cbind(exprs.pop.LRT.PC1, popQ=popQ$qvalues, popQsig=popQ$significant)
#write pval table. slow? Be sure to include cov name
write.table(exprs.pop.LRT.PC1.Q, file="lme4_qval_PC1_pop.txt", sep="\t")

###
#add to 1st PC1 table
PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") 
PC1qpop <- read.table("lme4_qval_PC1_pop.txt", header=T, sep="\t") 

#check numbers
test <- merge(PC1q[,1:10], PC1qdr)
write.table(test, file="lme4_qval_PC1.txt", sep="\t")