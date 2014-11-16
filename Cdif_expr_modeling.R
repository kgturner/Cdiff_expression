#Statistical modeling run on cluster (darjeeling) using maanova package
#11/11/2014

#Install/compile/add to path R version 3.1.1

#set up screen in shell
#screen

#start R session in shell
#R

#install necessary packages
source("http://bioconductor.org/biocLite.R")
#install packages
biocLite(c("limma", "oligo", "maanova", "pdInfoBuilder"))
#biocLite("limma")
#biocLite("maanova")
#biocLite("pdInfoBuilder")
biocLite("qvalue")

#load libraries
library(limma)
library(oligo)
library(maanova)
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

# Create a madata object #data is already normalized, background corrected, summarized. don't log transform.
madata <- read.madata(exprs.Cdifxys, exprs.design, log.trans=F)

###########Mixed model, pop and sampleID#####
# # Fit the model #there may be a warning message here

fit.full.mix <- fitmaanova(madata, formula=~Origin*Trt+ Tmpt + Pop+SampleID,covariate=~Tmpt, random=~Pop+SampleID, verbose=TRUE) #fit model for treatment, treatment set (timepoints), population type. everything else is 'random effect'. fixed effect model. we can also specify colums to test as random effects for a random effect model. 

#one permutation, can i test all these factors
test.full.mix.1=matest(madata, fit.full.mix, term=c("Origin", "Trt", "Tmpt"), n.perm=1, shuffle.method="sample", verbose=TRUE)

ftest.all = matest(madata, fit.full.mix, test.method=c(1,1),shuffle.method="sample", term=c("Origin*Trt", "Origin","Trt","Tmpt"), n.perm= 100)

ftest.adjP = adjPval(ftest.all, method = 'jsFDR')

summarytable(ftest.adjP)

pdf("fullmix_adjP.pdf")
volcano(ftest.adjP)
dev.off()


idx.full.mix <- volcano(ftest.all,method=c("unadj", "fdrperm"),highlight.flag=FALSE)

#clustering
cluster.kmean <- macluster(fit.full.mix, term="Origin",idx.gene=idx.full.mix$comparison1$idx.Fs,what="gene", method="kmean",kmean.ngroups=5, n.perm=100)

con.kmean <- consensus(cluster.kmean, 0.7)
con.kmean$groupname
cluster.hc <- macluster(fit.full.mix, term="Origin",idx.gene=idx.full.mix$comparison1$idx.Fs,what="sample", method="hc", n.perm=100)
con.hc <- consensus(cluster.hc)

# #test sample effect
 # test.samp <- matest(madata, fit.all, term="SampleID", n.perm=1)
 # #pval for sample
 # library(qvalue)
 # test.samp.p<- adjPval(test.samp, method=jsFDR)
 
 # #save summary table
# result = summarytable(test.all.1)

# #volcano plot
# pdf(file="test_all_1_volcano.pdf")
# idx.fit.all <- volcano(test.all.1,method=c("unadj", "fdrperm"),
# highlight.flag=FALSE)
# dev.off()
 
# #clustering
# cluster.kmean <- macluster(fit.all, ,term="SampleID",idx.gene=idx.fit.all$idx.all,what="gene", method="kmean",kmean.ngroups=5, n.perm=100)
# con.kmean <- consensus(cluster.kmean, 0.7)
# con.kmean$groupname
 
# cluster.hc <- macluster(fit.all, term="SampleID",idx.gene=idx.fit.all$idx.all,what="sample",method="hc", n.perm=100)
# con.hc <- consensus(cluster.hc) 
 
 
#########Mixed model, pop only###########################
# Fit the model #there may be a warning message here
#### maybe sampleID not necessary in model? since no technical replicates?
fit.mix <- fitmaanova(madata, formula=~Origin*Trt+ Tmpt + Pop,covariate=~Tmpt, random=~Pop, verbose=TRUE) 
#fit model for treatment, treatment set (timepoints), population type. everything else is 'random effect'. mixed effect model. we can also specify colums to test as random effects for a random effect model. 

#one permutation, can i test all these factors
test.mix.1=matest(madata, fit.mix, term=c("Origin", "Trt", "Tmpt"), n.perm=1, shuffle.method="sample", verbose=TRUE)

ftest.all = matest(madata, fit.mix, test.method=c(1,1),shuffle.method="sample", term=c("Origin","Trt","Tmpt"), n.perm= 100)

ftest.adjP = adjPval(ftest.all, method = 'jsFDR')

summarytable(ftest.adjP)

pdf("mix_adjP.pdf")
volcano(ftest.adjP)
dev.off()


idx.mix <- volcano(ftest.all,method=c("unadj", "fdrperm"),highlight.flag=FALSE)

#clustering
cluster.kmean <- macluster(fit.mix, term="Origin",idx.gene=idx.mix$comparison1$idx.Fs,what="gene", method="kmean",kmean.ngroups=5, n.perm=100)

con.kmean <- consensus(cluster.kmean, 0.7)
con.kmean$groupname
cluster.hc <- macluster(fit.mix, term="Origin",idx.gene=idx.mix$comparison1$idx.Fs,what="sample", method="hc", n.perm=100)
con.hc <- consensus(cluster.hc)

#########Fixed model (can do locally)###########################
# Fit the model #there may be a warning message here
fit.fix <- fitmaanova(madata, formula=~Origin*Trt+ Tmpt,covariate=~Tmpt, verbose=TRUE)

#fit model for treatment, treatment set (timepoints), population type. everything else is 'random effect'. 
#one permutation, can i test all these factors
test.fix.1=matest(madata, fit.fix, term=c("Origin", "Trt", "Tmpt"), n.perm=1, shuffle.method="sample", verbose=TRUE)

ftest.all = matest(madata, fit.fix, test.method=c(1,1),shuffle.method="sample", term=c("Origin","Trt","Tmpt"), n.perm= 100)
summarytable(ftest.all,outfile="fix_summarytable.csv")

ftest.adjP = adjPval(ftest.all, method = 'jsFDR')
summarytable(ftest.adjP, outfile="fix_adjP_summarytable.csv")

pdf("fix_adjP.pdf")
volcano(ftest.adjP)
dev.off()


idx.fix <- volcano(ftest.adjP,method=c("unadj", "fdrperm"),highlight.flag=FALSE)




#clustering
cluster.kmean <- macluster(fit.fix, term="Origin",idx.gene=idx.fix$comparison1$idx.Fs,what="gene", method="kmean",kmean.ngroups=5, n.perm=100)

con.kmean <- consensus(cluster.kmean, 0.7)
con.kmean$groupname
cluster.hc <- macluster(fit.fix, term="Origin",idx.gene=idx.fix$comparison1$idx.Fs,what="sample", method="hc", n.perm=100)
con.hc <- consensus(cluster.hc)





#############Evans notes#########
# #test interaction of factors
# #population_type_treatment_set means & standard errors
# population_type_treatment_set.means = t(apply(madata$data, 1, tapply, design$population_type_treatment_set, mean)) 
# colnames(population_type_treatment_set.means)=paste("Mean", colnames(population_type_treatment_set.means), sep=":")
# population_type_treatment_set.SEs = t(apply(madata$data, 1, tapply, design$population_type_treatment_set, function(x) sqrt(var(x)/length(x))))
# colnames(population_type_treatment_set.SEs)=paste("SE", colnames(population_type_treatment_set.SEs), sep=":")

# #differences in means
# #don't worry about this for now


# #find interaction genes
# test.int=matest(madata, fit.fix, term=c("population_type_treatment_set"), n.perm=1000, shuffle.method="sample", verbose=TRUE)  #leaving off here Jan 9th afternoon
# library(qvalue)
# test.int=adjPval(test.int, method='jsFDR') #pval adjust

# #and non-interaction genes

############lme4 modeling?#####
#use design and exprs.Cdifxys
#exprs.Cdifxys needs to be rotated, and desing cols added? but what is the response variable?


modelfull <- lmer()
model1<-lmer(LfLgth1  ~ Origin* Latitude +(1|PopID/MomFam), family=gaussian,data=modeldata) #Includes maternal family variance

