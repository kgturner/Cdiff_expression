#STEPS ALL TAKE PLACE ON MRFOX NOW

#steps
#1. identify bad areas (using NimbleScan, as recorded in microarray_QC.R), replace their values with NAs in the .xys files (using bad_features.kay_microarray.pl).
#2. load data into R, impute values using impute.knn (this R file)
#3. export imputed values, replace NAs in .xys files with imputed values (using knnreplace.kay_microarray.pl)
#4. load modified .xys files into R using oligo, prune bad samples, run rma (this R file)
#5. ready for analysis (as in Kay_analysis.txt or Kay_Final_analysis.txt - actually R code)

#specify bioconductor source
source("http://bioconductor.org/biocLite.R")
#install packages
biocLite(c("affy"))

#loaded libraries: oligo, limma
#library(limma)
#library(affy)
library(oligo)
library(impute)
#library(genefilter)

############################################################
############################################################
#construct annotation library
library(pdInfoBuilder)
baseDir <- "/home/morien/microarray/sunflower/tif_aligned/corrected_xys"
xys <- list.files(baseDir, pattern = ".xys", full.names = TRUE)[1]
baseDir <- "/home/morien/microarray/sunflower/design_files"
ndf <- list.files(baseDir, pattern = ".ndf", full.names = TRUE)
seed <- new("NgsExpressionPDInfoPkgSeed", ndfFile = ndf, xysFile = xys, biocViews = "AnnotationData", organism = "sunflower", species = "Helianthus annuus")
makePdInfoPackage(seed, destDir = baseDir)
#now install package. there is no need to load it in after. R will recognize that it's installed once you try to read the xys files
install.packages("/home/morien/microarray/sunflower/design_files/pd.110210.hannuus.lz.exp/", repos=NULL, type="source")
############################################################
############################################################

#read in the data #working on mrfox
xys.files<-list.xysfiles("/home/morien/microarray/AMF/Helianthus_amf_TIF_aligned/corrected_xys", full.names=TRUE)
xys.files2<-list.xysfiles("/home/morien/microarray/AMF/Helianthus_amf_TIF_aligned/19933_20090409/corrected_xys", full.names=TRUE)
xys.files=c(xys.files, xys.files2)
AMF<-read.xysfiles(xys.files)
xys.files<-list.xysfiles("/home/morien/microarray/NIL/tif_aligned/corrected_xys", full.names=TRUE)
NIL<-read.xysfiles(xys.files)
xys.files<-list.xysfiles("/home/morien/microarray/sunflower/tif_aligned/corrected_xys", full.names=TRUE)
sunflower<-read.xysfiles(xys.files)
xys.files<-list.xysfiles("/home/morien/microarray/star_thistle/tif_aligned/corrected_xys", full.names=TRUE)
star_thistle<-read.xysfiles(xys.files)


#remove random probes	#identifies probes which only have NA values across every array
data<-exprs(NIL)
data2<-data[apply(data,1,function(x)any(!is.na(x))),]
data3<-data[apply(data,1,function(x)any(is.na(x))),]
data2<-data[rowSums(is.na(data)) < 129, ,drop=FALSE] #number here should match the number of samples (columns) in the data frame

#k nearest neighbor imputation
knnout<-impute.knn(data2, k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=36243606) #imputing with default values
allkn<-knnout$data

#print out and replace original values with imputed values using kay's script (knnreplace.kay_mircroarray)
write.table(allkn, file="microarray/NIL/knn_replaced.txt", quote=F, sep="\t")

#look at raw data
#hist(data2, col = darkColors(132), lty = 1, xlim = c(6, 17))
#boxplot(data2[, 1:132], col = darkColors(132), names = 1:132)

#read in KNN imputed XYS files
xys.files<-list.xysfiles("/home/morien/microarray/AMF/Helianthus_amf_TIF_aligned/knn_xys", full.names=TRUE)
xys.files2<-list.xysfiles("/home/morien/microarray/AMF/Helianthus_amf_TIF_aligned/19933_20090409/knn_xys", full.names=TRUE)
xys.files=c(xys.files, xys.files2)
AMF<-read.xysfiles(xys.files)
xys.files<-list.xysfiles("/home/morien/microarray/NIL/tif_aligned/knn_xys", full.names=TRUE)
NIL<-read.xysfiles(xys.files)
xys.files<-list.xysfiles("/home/morien/microarray/sunflower/tif_aligned/knn_xys", full.names=TRUE)
sunflower<-read.xysfiles(xys.files)
xys.files<-list.xysfiles("/home/morien/microarray/star_thistle/tif_aligned/knn_xys", full.names=TRUE)
star_thistle<-read.xysfiles(xys.files)

#RMA - background correction, normalization, summarization of probes
ppData.sunflower=rma(sunflower)
ppData.NIL=rma(NIL)
ppData.AMF=rma(AMF)
ppData.star_thistle=rma(star_thistle)

#RMA is complete, now extract expression data
exprs.sunflower=exprs(ppData.sunflower)
exprs.NIL=exprs(ppData.NIL)
exprs.AMF=exprs(ppData.AMF)
exprs.star_thistle=exprs(ppData.star_thistle)


#now remove genes that are unexpressed
#find random probe average expression and standard deviation
a=read.delim("microarray/random_probe_positions.txt", sep="\t", header=F) #one array should be enough for a test. we can use all the arrays if need be
b=log10(a[,10]) #column 10 of this file has the expression values for the random probes. we are comparing to log 10 transformed values, so transform random exp values too
> mean(b)
[1] 2.60122 #average random probe exp
> sd(b)
[1] 0.05576136 #random probe standard deviation

#remove 'unexpressed' probes

#read in raw data with probe info

#extracting probe expression values from sunflower probe info files #DONE IN SHELL NOT R
awk '{a[FNR]=a[FNR] FS $8;t=(FNR>T)?FNR:t}END {for (i=1;i<=t;i++) print a[i]}' 470975A01_532.probe 471139A08_532.probe 471148A03_532.probe 471186A10_532.probe 471189A05_532.probe 471198A12_532.probe 471225A07_532.probe 470975A02_532.probe 471139A09_532.probe 471148A04_532.probe 471186A11_532.probe 471189A06_532.probe 471208A01_532.probe 471225A08_532.probe 470975A03_532.probe 471139A10_532.probe 471148A05_532.probe 471186A12_532.probe 471189A07_532.probe 471208A02_532.probe 471225A09_532.probe 470975A04_532.probe 471139A11_532.probe 471148A06_532.probe 471187A01_532.probe 471189A08_532.probe 471208A03_532.probe 471225A10_532.probe 470975A05_532.probe 471139A12_532.probe 471148A07_532.probe 471187A02_532.probe 471189A09_532.probe 471208A04_532.probe 471225A11_532.probe 470975A06_532.probe 471145A01_532.probe 471148A08_532.probe 471187A03_532.probe 471189A10_532.probe 471208A05_532.probe 471225A12_532.probe 470975A07_532.probe 471145A02_532.probe 471148A09_532.probe 471187A04_532.probe 471189A11_532.probe 471208A06_532.probe 471312A01_532.probe 470975A08_532.probe 471145A03_532.probe 471148A10_532.probe 471187A05_532.probe 471189A12_532.probe 471208A07_532.probe 471312A02_532.probe 470975A09_532.probe 471145A04_532.probe 471148A11_532.probe 471187A06_532.probe 471198A01_532.probe 471208A08_532.probe 471312A03_532.probe 470975A10_532.probe 471145A05_532.probe 471148A12_532.probe 471187A07_532.probe 471198A02_532.probe 471208A09_532.probe 471312A04_532.probe 470975A11_532.probe 471145A06_532.probe 471186A01_532.probe 471187A08_532.probe 471198A03_532.probe 471208A10_532.probe 471312A05_532.probe 470975A12_532.probe 471145A07_532.probe 471186A02_532.probe 471187A09_532.probe 471198A04_532.probe 471208A11_532.probe 471312A06_532.probe 471139A01_532.probe 471145A08_532.probe 471186A03_532.probe 471187A10_532.probe 471198A05_532.probe 471208A12_532.probe 471312A07_532.probe 471139A02_532.probe 471145A09_532.probe 471186A04_532.probe 471187A11_532.probe 471198A06_532.probe 471225A01_532.probe 471312A08_532.probe 471139A03_532.probe 471145A10_532.probe 471186A05_532.probe 471187A12_532.probe 471198A07_532.probe 471225A02_532.probe 471312A09_532.probe 471139A04_532.probe 471145A11_532.probe 471186A06_532.probe 471189A01_532.probe 471198A08_532.probe 471225A03_532.probe 471312A10_532.probe 471139A05_532.probe 471145A12_532.probe 471186A07_532.probe 471189A02_532.probe 471198A09_532.probe 471225A04_532.probe 471312A11_532.probe 471139A06_532.probe 471148A01_532.probe 471186A08_532.probe 471189A03_532.probe 471198A10_532.probe 471225A05_532.probe 471312A12_532.probe 471139A07_532.probe 471148A02_532.probe 471186A09_532.probe 471189A04_532.probe 471198A11_532.probe 471225A06_532.probe > datamatrix.txt

sunprobes=read.delim("microarray/sunflower/tif_aligned/probe_reports/datamatrix.txt", sep=" ", header=T, row.names=NULL)

#above output doesn't contain probe names, get them from a single probe file
a=read.delim("microarray/sunflower/tif_aligned/probe_reports/470975A01_532.probe", sep="\t", header=T, row.names=NULL, skip=1) #duplicate probe ID names means we have to make them unique somehow
b=do.call(paste, c(a[c("PROBE_ID", "MATCH_INDEX")], sep = "_"))	#concatenating the match index with the probe ID produces a human readable unique probe ID
row.names(sunprobes)=b #assign unique row names to data matrix

#remove extra column introduced by awk
sunprobes=sunprobes[,-1]

#REMOVE NON-RANDOM CONTROL PROBES FOR "BETTER" NORMALIZATION AND RMA
index1 <- with(sunprobes, grepl("CPK", row.names(sunprobes))) #which probes are CPK probes
index2 <- with(sunprobes, grepl("AMPFS", row.names(sunprobes))) #which probes are other non-random control probes
index3 <- with(sunprobes, grepl("DELFS", row.names(sunprobes))) #which probes are other non-random control probes
index4 <- with(sunprobes, grepl("XENOTRACK", row.names(sunprobes))) #which probes are other non-random control probes
index5 <- with(sunprobes, grepl("EMPTY", row.names(sunprobes))) #which probes are other non-random control probes
index6 <- with(sunprobes, grepl("Syn[A-Z]_", row.names(sunprobes))) #which probes are other non-random control probes
a=which(index1 == "FALSE") #I WANT TO REMOVE THESE
b=which(index2 == "FALSE") #I WANT TO REMOVE THESE
c=which(index3 == "FALSE") #I WANT TO REMOVE THESE
d=which(index4 == "FALSE") #I WANT TO REMOVE THESE
e=which(index5 == "FALSE") #I WANT TO REMOVE THESE
f=which(index6 == "FALSE") #I WANT TO REMOVE THESE
g=Reduce(intersect, list(a,b,c,d,e,f)) #only want those that are false in all cases
test=sunprobes[g,]


#RMA with limma this time. no using oligo without .xys inputs. 
#background correction
bgc=backgroundCorrect(test, method="normexp", normexp.method="rma", verbose=T)
#quantile normalization
QNsunprobes=normalizeBetweenArrays(bgc, method="quantile")

#just gets the random probes on their own
#test=subset(QNsunprobes, grepl("RANDOM", row.names(QNsunprobes))) #just gets the random probes

#define random probes by index
QNsunprobes=as.data.frame(QNsunprobes)
index1 <- with(QNsunprobes, grepl("RANDOM", row.names(QNsunprobes))) #which probes are RANDOM probes
a=which(index1 == "FALSE") #these are the non-random probes
b=which(index1 == "TRUE") #these are the random probes
rand=QNsunprobes[b,]
data=QNsunprobes[a,]

#unexpressed <= 2 * (Standard deviation of random probes) + (mean random probes)
> (0.05576136 * 2) + 2.60122
[1] 2.712743	#remove genes with expression below this level

min(exprs.sunflower)
[1] 3.959489	#no genes to remove

#check untransformed data (diff between random and actual probes)

#using MAANOVA package to analyze data w a mixed model
#building this from MAANOVA tutorial

# Read hybridization design
design=read.table("microarray/sunflower/design_files/experimental_design.txt", header=T, sep="\t") #sunflower expt design is located in Sunflower_HX12_microarrayFinalInfo.xls file, transferred to microarray/sunflower/design_files/experimental_design.txt and modified to remove bad samples and add additional metadata

# show design
print(design)

#load maanova library
library(maanova)

# Create a madata object #data is already normalized, background corrected, summarized. don't log transform.
madata=read.madata(exprs.sunflower, design, log.trans=F)

#########################################################################################################################
#design:
#
#treatment_set (jas1 jas2 dr1 dr2 con)
#population_type (weedy wild dom)
#treatment_set * population_type
#population (random)

#1. test the interaction first
#2. then test the main effects treatment_set and population type on the remaining genes that do not have an interaction
#3. for the interaction genes compare weedy vs wild, weedy vs dom, and dom vs wild for each treatment_set
#4. for the significant population_type genes compare weedy vs wild, weedy vs dom, and dom vs wild across all treatments
#5. for the significant treatment genes do all the pairwise comparisons between the treatments
#########################################################################################################################

# Fit the model #there may be a warning message here
fit.fix=fitmaanova(madata, formula=~treatment_set+population_type+population_type_treatment_set+population, random=~population) #fit model for treatment, treatment set (timepoints), population type. everything else is 'random effect'. fixed effect model. we can also specify colums to test as random effects for a random effect model. 

#one permutation, can i test all these factors
test.one=matest(madata, fit.fix, term=c("treatment_set", "population_type", "population_type_treatment_set"), n.perm=1, shuffle.method="sample", verbose=TRUE)

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
#check kay's file for these


#separate non-interaction genes

#test those genes for treatment_set & population_type
test.TS.PT=matest(madata, fit.fix, term=c("treatment_set", "population_type"), n.perm=1, shuffle.method="sample", verbose=TRUE)


#1000 permutations
test.full=matest(madata, fit.fix, term=c("treatment_set", "population_type", "population_type_treatment_set"), n.perm=1000, shuffle.method="sample",verbose=TRUE)  


#find interactions between experimental groups
#treatment_set
treatment_set.means = t(apply(madata$data, 1, tapply, design$treatment_set, mean))
colnames(treatment_set.means)=paste("Mean", colnames(treatment_set.means), sep=":")
treatment_set.SEs = t(apply(madata$data, 1, tapply, design$treatment_set, function(x) sqrt(var(x)/length(x))))
colnames(treatment_set.SEs)=paste("SE", colnames(treatment_set.SEs), sep=":")

#population_type
population_type.means = t(apply(madata$data, 1, tapply, design$population_type, mean)) 
colnames(population_type.means)=paste("Mean", colnames(population_type.means), sep=":")
population_type.SEs = t(apply(madata$data, 1, tapply, design$population_type, function(x) sqrt(var(x)/length(x))))
colnames(population_type.SEs)=paste("SE", colnames(population_type.SEs), sep=":")




#pval correction
library(qvalue)
test.full.adj=adjPval(test.full, method='jsFDR')

#calculate differences in means between different groups

#population_type
domvsweedy=matrix(population_type.means[,1]-population_type.means[,2])
domvswild=matrix(population_type.means[,1]-population_type.means[,3])
weedyvswild=matrix(population_type.means[,2]-population_type.means[,3])
domvsweedypluswild=matrix(population_type.means[,1]-((population_type.means[,2]+population_type.means[,3])/2))

#treatment
ctlvsjas=matrix(treatment.means[,1]-treatment.means[,2])
ctlvsdr=matrix(treatment.means[,1]-treatment.means[,3])
jasvsdr=matrix(treatment.means[,2]-treatment.means[,3])
ctlvsjasplusdr=matrix(treatment.means[,1]-((treatment.means[,2]+treatment.means[,3])/2))

#treatment_set
jas1vsjas2=matrix(treatment_set.means[,4]-treatment_set.means[,5])
dr1vsdr2=matrix(treatment_set.means[,2]-treatment_set.means[,3])

#population
#not doing this right now. too many comparisons

#bind differences into matrix
meandiffs=cbind(domvsweedy, domvswild, weedyvswild, domvsweedypluswild, ctlvsjas, ctlvsdr, jasvsdr, ctlvsjasplusdr, jas1vsjas2, dr1vsdr2)

colnames(meandiffs)=c("domvsweedy", "domvswild", "weedyvswild", "domvsweedypluswild", "ctlvsjas", "ctlvsdr", "jasvsdr", "ctlvsjasplusdr", "jas1vsjas2", "dr1vsdr2")

row.names(meandiffs)=row.names(exprs.sunflower)

#make contrast matrices
cmat.population_type=makeContrasts(
domestic-weedy, domestic-wild, weedy-wild, domestic - (weedy + wild),
levels=c("weedy", "domestic", "wild"))
cmat.population_type=t(cmat.population_type)

#cmat.treatment=makeContrasts(
#CTL - JAS, CTL - DR, JAS-DR, CTL - (JAS + DR),
#levels=c("CTL", "JAS", "DR"))
#cmat.treatment=t(cmat.treatment)

cmat.treatment_set=makeContrasts(
#CTL - JAS1, CTL - JAS2, CTL - DR1, CTL - DR2, JAS1-JAS2, DR1-DR2,
JAS1 - DR2, JAS2 - DR1,
levels=c("CTL", "JAS1", "JAS2", "DR1", "DR2"))
cmat.treatment_set=t(cmat.treatment_set)

#run model with contrast matrices
fit.fix=fitmaanova(madata, formula=~treatment+treatment_set+population_type+sample+population, random=~sample+population)
test.full.pop_type=matest(madata, fit.fix, term="population_type", Contrast=cmat.population_type, n.perm=1, shuffle.method="sample",verbose=TRUE)	#the number 4 test is not estimable (?!)
#test.full.treatment=matest(madata, fit.fix, term="treatment", Contrast=cmat.treatment, n.perm=1, shuffle.method="sample",verbose=TRUE)			#none of these are estimable ?????
test.full.treatment_set=matest(madata, fit.fix, term="treatment_set", Contrast=cmat.treatment_set, n.perm=1, shuffle.method="sample",verbose=TRUE)	#this works fine (!?)
#no idea what is going on here



#boxplot of diffs, ordered by maximum value
absmeandiffs=abs(meandiffs)
test1=apply(absmeandiffs, 2, q3)

boxplot(abs(meandiffs[,order(test1)]), las=2)


###################### #function to transform log diffs into fold change diffs
logdiff2FC=function(x, base=2) {
  # Transforms x from a log difference (logdiff) to a fold change (FC) scale
  # logdiff is the expression difference in log scale of base = 'base'
  # FC is the ratio between two expression values in the original scale
  # The absolute value of the returned FC is always greater than or equal to 1
  # The sign is kept to indicate direction of change
  # When the the log difference is 0, the FC is always 1.
  signs = sign(x)
  out   = base^abs(x)
  if(any(x==0)) {
    out[x==0] = 1
  }
  return(out*signs)
} # logdiff2FC
######################


# Transform log differences to FC scale
FC = apply(meandiffs, 2, logdiff2FC)



#make colnames labels more specific
colnames(FC)=paste("FC", colnames(FC), sep=":")
colnames(meandiffs)=paste("MD", colnames(meandiffs), sep=":")

out = data.frame(treatment.means, treatment.SEs, treatment_set.means, treatment_set.SEs, population_type.means, population_type.SEs, F1_val=test.full.adj$F1$Fobs, Fs_val=test.full.adj$Fs$Fobs, P_tab=test.full.adj$Fs$Ptab, FDR_tab=test.full.adj$Fs$adjPtab, meandiffs, FC)

#which genes are significant after FDR correction






######################
q3 <- function (x) { #function to compute the third quartile of a vector
    x<- sort(x)
    n <- length(x)
    quot <- 3*(n+1) %/% 4
    rem <- 3*(n+1) %% 4
    if(rem == 0 | rem == 1) {
        return (x[quot])
    }
    else if(rem ==2 | rem == 3) {
        return (0.5 * (x[quot] + x[quot+1]))
    }
}
######################


#make the contrast matrix
#could do by hand... this will take ages
cmat = rbind( Geno     =  c( 1,   1,  -1,  -1 )*.5,
              Trt      =  c( 1,  -1,   1,  -1 )*.5,

#limma has a makeContrasts function #remember that reverse comparisons are not equivalent

cmat.population_type=makeContrasts(
domestic-weedy, domestic-wild, weedy-wild, domestic-(weedy+wild), weedy-(domestic+wild), wild-(weedy+domestic),
levels=c("weedy", "domestic", "wild"))

cmat.treatment=makeContrasts(
CTL-JAS, CTL-DR, JAS-DR, CTL-(JAS+DR), 
levels=c("CTL", "JAS", "DR"))

cmat.treatment_set=makeContrasts(
JAS1-JAS2, DR1-DR2,
levels=c("CTL", "JAS1", "JAS2", "DR1", "DR2"))
))

cmat.population=makeContrasts( #mostly the 1 to 1 population effects will not be particularly interesting. i'll code these up if we think there are any specific comparisons we should make

levels=c())

cmat.treatment_group = makeContrasts(
#weedy vs wild vs domestic
a1-b1, a1-c1, b1-c1, a2-b2, a2-c2, b2-c2, a3-b3, a3-c3, b3-c3, a1-(b1+c1)/2, a2-(b2+c2)/2, a3-(b3+c3)/2,

#within-population type comparisons
a1-a2, a1-a3, a2-a3, a1-(a2+a3)/2, b1-b2, b1-b3, b2-b3, b1-(b2+b3)/2, c1-c2, c1-c3, c2-c3, c1-(c2+c3)/2,

levels=c("a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3"))

#levels	
#a1=domestic_ctl
#a2=domestic_jas
#a3=domestic_dr
#b1=weedy_ctl
#b2=weedy_jas
#b3=weedy_dr
#c1=wild_ctl
#c2=wild_jas
#c3=wild_dr

dim(cmat.treatment_group)
[1] 9 60 #nine groups by 60 comparisons

# Log differences for main effects and the interaction
effect = apply(Means, 1, function(x) sum(x*cmat.treatment_group[,n])) #where means is the object that shows the mean expression for each treatment group and n is the row of the cmat matrix which delineates the effect we are interested in
effect1 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,1])) #a1-b1 #domestic vs weedy controls
effect2 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,2])) #"a1 - c1" #domestic vs wild controls
effect3 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,3])) #"b1 - c1" #weedy vs wild controls
effect4 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,4])) #"a2 - b2" #domestic vs weedy JAS
effect5 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,5])) #"a2 - c2" #domestic vs wild JAS
effect6 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,6])) #"b2 - c2" #weedy vs wild JAS
effect7 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,7])) #"a3 - b3" #domestic vs weedy DR
effect8 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,8])) #"a3 - c3" #domestic vs wild DR
effect9 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,9])) #"b3 - c3" #weedy vs wild DR
effect10 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,10])) #"a1 - (b1 + c1)/2" #domestic vs weedy+wild controls
effect11 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,11])) #"a2 - (b2 + c2)/2" #domestic vs weedy+wild JAS
effect12 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,12])) #"a3 - (b3 + c3)/2" #domestic vs weedy+wild DR
effect13 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,13])) #"b1 - a1" #weedy vs domstic controls
effect14 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,14])) #"c1 - a1" #wild vs domestic controls
effect15 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,15])) #"c1 - b1" #wild vs weedy controls
effect16 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,16])) #"b2 - a2" #weedy vs domestic JAS
effect17 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,17])) #"c2 - a2" #wild vs domestic JAS
effect18 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,18])) #"c2 - b2" #wild vs weedy JAS
effect19 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,19])) #"b3 - a3" #weedy vs domestic DR
effect20 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,20])) #"c3 - a3" #wild vs domestic DR
effect21 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,21])) #"c3 - b3" #wild vs weedy DR
effect22 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,22])) #"((b1 + c1)/2) - a1" #weedy&wild vs domestic controls
effect23 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,23])) #"((b2 + c2)/2) - a2" #weedy&wild vs domestic JAS
effect24 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,24])) #"((b3 + c3)/2) - a3" #weedy&wild vs domestic DR
effect25 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,25])) #"a1 - a2" #domestic controls vs JAS
effect26 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,26])) #"a1 - a3" #domestic controls vs DR
effect27 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,27])) #"a2 - a3" #domestic JAS vs DR
effect28 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,28])) #"a1 - (a2 + a3)/2" #domestic controls vs treatments
effect29 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,29])) #"b1 - b2" #weedy controls vs JAS
effect30 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,30])) #"b1 - b3" #weedy controls vs DR
effect31 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,31])) #"b2 - b3" #weedy JAS vs DR
effect32 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,32])) #"b1 - (b2 + b3)/2" #weedy controls vs treatments
effect33 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,33])) #"c1 - c2" #wild controls vs JAS
effect34 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,34])) #"c1 - c3" #wild controls vs DR
effect35 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,35])) #"c2 - c3" #wild JAS vs DR
effect36 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,36])) #"c1 - (c2 + c3)/2" #wild controls vs treatments
effect37 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,37])) #"a2 - a1" #domestic JAS vs controls
effect38 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,38])) #"a3 - a1" #domestic DR vs controls
effect39 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,39])) #"a3 - a2" #domestic DR vs JAS
effect40 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,40])) #"((a2 + a3)/2) - a1" #domestic treatments vs controls
effect41 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,41])) #"b2 - b1" #weedy JAS vs controls
effect42 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,42])) #"b3 - b1" #weedy DR vs controls
effect43 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,43])) #"b3 - b2" #weedy DR vs JAS
effect44 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,44])) #"((b2 + b3)/2) - b1" #weedy treatments vs controls
effect45 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,45])) #"c2 - c1" #wild JAS vs controls
effect46 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,46])) #"c3 - c1" #wild DR vs controls
effect47 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,47])) #"c3 - c2" #wild DR vs JAS
effect48 = apply(means, 1, function(x) sum(x*cmat.treatment_group[,48])) #"((c2 + c3)/2) - c1" #wild JAS&DR vs controls

# Bind columns into a matrix #excluded: effect52,
logDiffs=cbind(effect1, effect2, effect3, effect4, effect5, effect6, effect7, effect8, effect9, effect10, effect11, effect12, effect13, effect14, effect15, effect16, effect17, effect18, effect19, effect20, effect21, effect22, effect23, effect24, effect25, effect26, effect27, effect28, effect29, effect30, effect31, effect32, effect33, effect34, effect35, effect36, effect37, effect38, effect39, effect40, effect41, effect42, effect43, effect44, effect45, effect46, effect47, effect48)

# Transform log differences to FC scale
FC = apply(logDiffs, 2, logdiff2FC)

cmat.treatment_group=t(cmat.treatment_group)
# Test test each contrasts using 300 permutations of sample labels #transpose cmat, the limma function to make it has the columns and rows reversed compared to maanova
test.cmat.treatment_group=matest(madata, fit.fix, term="treatment_group", Contrast=cmat.treatment_group, n.perm=300, test.type = "ttest",
                 shuffle.method="sample", verbose=TRUE)

# Contrasts names are not kept in the matrix of permutation results, so lets 
# copy them from the matrix of tabular p-values
colnames(test.cmat.treatment_group$Fs$Pvalperm) = colnames(test.cmat.treatment_group$Fs$Ptab)

# Multiple comparison control (FDR transformation, see ?adjPval)
test.cmat.treatment_group=adjPval(test.cmat.treatment_group, method="adaptive")




#eliminate reverse contrasts

#

#test for interaction effects

#test the contrasts for groups where we expect an interaction effect

#where there is a significant treatment effect, test for sig diff genes















############################################################
############################################################
#LIMMA RMA method (alternative to oligo RMA). only use if necessary.
#background correction
bgc=backgroundCorrect(allkn, method="normexp", normexp.method="rma", verbose=T)
#quantile normalization
QNsunflower=normalizeBetweenArrays(bgc, method="quantile")

############################################################
############################################################
#kay's code
#KEEP random probes
#knn - this imputes data for the bad features
knnout<-impute.knn(data2, k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=36243606)
allkn<-knnout$data
subknall<-allkn[,c("463790A02_532.xys","463778A10_532.xys", "463785A01_532.xys",  "463790A12_532.xys", "463884A01_532.xys", "463884A03_532.xys", "463884A01_532.xys", "463968A03_532.xys", "463968A05_532.xys", "463968A10_532.xys")]
write.table(subknall, file="knn")

use knnreplace.pl to substitute knn values for bad values

#rawdata
hist(ragc, col = darkColors(59), lty = 1, xlim = c(6, 30))
boxplot(ragc[, 1:143], col = darkColors(143), names = 1:143)

#knn corrected raw data
#preprocessing includes: background, normalizing, expression - can do in one step or separately 

#preprocessing in one step
ppData<-rma(ragc)

#create boxplots
boxplot(ppData[, 1:143], col = darkColors(143), names = 1:143)
hist(ppData, col = darkColors(59), lty = 1, xlim = c(6, 30))

#background
bData<-backgroundCorrect(ragc, method="rma")
#normalization alone
bnData<-normalize(exprs(bData), method="quantile")

#summarize using probe IDs:
bnsData<-summarize(noNA, probes=infor[,1], method="medianpolish")


#how to screen for genes with low expression across the board
 
t<-(y[apply(y,1,function(x)any(x>9.453454))])
tm<-matrix(t, ,45062)
td<-data.frame(tm)


#I removed probes with expression levels within 2 SD of the random probe normalized intensity values
#cutoff is 95% 4 out of 143 are above cutoff #double checked if random probes right and if number without random correspond to oligo ram->everything is right
#cut-off= SD*2+mean
0.5195944*2+ 8.414265  
#=9.453454

col<-ncol(bnsData)

number<-rowSums(matrix(as.integer(bnsData>9.453454), ,col))
fexprs<-bnsData[number>7, ]
 
dim(bnsData)
[1] 45063   143
> dim(fexprs)
[1] 33464   143
