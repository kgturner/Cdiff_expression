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
biocLite(c("limma", "oligo", "lme4", "pdInfoBuilder", "qvalue"))
#biocLite("limma")
#biocLite("maanova")
#biocLite("pdInfoBuilder")
# biocLite("qvalue")

#load libraries
library(limma)
library(oligo)
library(lme4)
library(pdInfoBuilder)
library(qvalue)

#now install custom package for this experiment. there is no need to load it in after. R will recognize that it's installed once you try to read the xys files
install.packages("~/Centaurea_diffusa_expression/pd.110405.cdiffusa.lz.exp/", repos=NULL, type="source")

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
write.table(exprs.clim, file="experimentaldesign.txt")

exprs.design <- read.table("C:/Users/Kat/Documents/GitHub/Cdiff_expression/experimentaldesign.txt", header=T, sep="\t") 

#add design info to data frame
exprs.df <- merge(exprs.design,exprs.df, all.y=TRUE)
#write data+design file. quite slow.
write.table(exprs.df, file="Cdifexprs_lme4dat.txt")
#read data table, probably slow
exprs.df <- read.table("C:/Users/Kat/Documents/GitHub/Cdiff_expression/Cdifexprs_lme4dat.txt", header=T, sep="\t") 

####loop writing####
# #for loop?
# for (i in colnames(exprs.df[,14:61037])){
#   model1 <- lmer(i ~ Origin*Trt+Latitude+ (Tmpt|SampleID)+(Origin|Pop), family=gaussian, data=exprs.df)
# }

#func?
exprs.LR<- function(trait,df,cov, family=gaussian){
  modeldata<-df[!is.na(df[[trait]]),]
  
  #browser()
  
  model1<-lmer(modeldata[[trait]]  ~ Origin*Trt+cov+ (Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)
  model2<-lmer(modeldata[[trait]]  ~ Origin+Trt+ cov + (Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)
  model3<-lmer(modeldata[[trait]]  ~ Origin+Trt + (Tmpt|PopTrtPool)+(1|Pop), family,data=modeldata)
  model4 <- lmer(modeldata[[trait]]  ~ Trt+ cov+(Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=modeldata)

  a1 <- anova(model2,model1) # is interaction sig?
  a2 <- anova(model3,model2) # is covariate sig?
  a3 <- anova(model4, model2) #is origin sig?
  
  pval <- as.data.frame(cbind(Contig=trait,intLRT=a1[[7]][2],covLRT=a2[[7]][2],originLRT=a3[[7]][2]))
  
  return(pval)
}
exprs.LRT <- lapply(names(exprs.df)[15:61038],function(n) exprs.LR(n,df=exprs.df,cov="Latitude"))#apply func to all things in list
row.names(exprs.LRT) <- names(exprs.df)[15:61038]

#testing
test <- exprs.df[c(1:6), c(1:20)]
# test$PopTrtPool <- paste0(test$Pop, "_",test$TrtPool)
model1<-lmer(Contig1  ~ Origin*Trt+Latitude+(Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=test)
summary(model1)
model2<-lmer(Contig1  ~ Origin+Trt+Latitude+(Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=test)
a1 <- anova(model2,model1)
model3<-lmer(Contig1  ~ Origin+Trt+ (Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=test)
a2 <- anova(model3,model2)
model4 <- lmer(Contig1  ~ Trt+ Latitude+(Tmpt|PopTrtPool)+(1|Pop), family=gaussian,data=test)
a3 <- anova(model4, model2)
#get out pval of LRT
a2[[7]][2]

#paste to new dataframe
#cols: contig, intLRT, covLRT, originLRT
#aka: trait, a1[[7]][2], ditto a2, ditto a3
trait <- "Contig1"
pval <- as.data.frame(cbind(Contig=trait,intLRT=a1[[7]][2],covLRT=a2[[7]][2],originLRT=a3[[7]][2]))

test <- exprs.df[c(1:6), c(1:20)]
test.LRT <- lapply(names(test)[15:20],function(n) exprs.LR(n,df=test,cov="Latitude"))#apply func to all things in list
