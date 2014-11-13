#Statistical modeling run on cluster (darjeeling)
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

#load libraries
library(limma)
library(oligo)
library(maanova)
library(pdInfoBuilder)

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

# Fit the model #there may be a warning message here
fit.fix <- fitmaanova(madata, formula=~Origin*Trt+ Tmpt + Pop+SampleID,covariate=~Tmpt, random=~Pop+SampleID, verbose=TRUE) #fit model for treatment, treatment set (timepoints), population type. everything else is 'random effect'. fixed effect model. we can also specify colums to test as random effects for a random effect model. 

#one permutation, can i test all these factors
test.one=matest(madata, fit.fix, term=c("Origin", "Trt", "Tmpt"), n.perm=1, shuffle.method="sample", verbose=TRUE)





#### maybe sampleID not necessary in fit.fix?









#test interaction of factors
#population_type_treatment_set means & standard errors
population_type_treatment_set.means = t(apply(madata$data, 1, tapply, design$population_type_treatment_set, mean)) 
colnames(population_type_treatment_set.means)=paste("Mean", colnames(population_type_treatment_set.means), sep=":")
population_type_treatment_set.SEs = t(apply(madata$data, 1, tapply, design$population_type_treatment_set, function(x) sqrt(var(x)/length(x))))
colnames(population_type_treatment_set.SEs)=paste("SE", colnames(population_type_treatment_set.SEs), sep=":")

#differences in means
#don't worry about this for now


#find interaction genes
test.int=matest(madata, fit.fix, term=c("population_type_treatment_set"), n.perm=1000, shuffle.method="sample", verbose=TRUE)  #leaving off here Jan 9th afternoon
library(qvalue)
test.int=adjPval(test.int, method='jsFDR') #pval adjust

#and non-interaction genes
