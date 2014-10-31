#Cdiffusa drought arrays preprocessing
#10/21/2014

#Preprocessing and analysis overview
#1. Grid .tif files and produce .xys file for each .tif and generate Experimental metrics report and Probe reports, using NimbleScan 2.6
#2. identify bad areas (coordinates from NimbleScan, recorded in Cdiff_expression_arrayQC.R, file names and coordinates for further processing in Cdiffexpression_bad_features.txt), 
#3. replace values of bad features with NAs in the .xys files (using bad_features.kay_microarray.pl; or batch call removebadfeatures.pl).
#4. Discard/move original .xys files with bad features. If more than one bad feature removed, retain most highly modified version of file (i.e., largest numbef of "_bfr" endings, then change file ending to .xys).
#5. load data into R, impute values using impute.knn (this R file)
#6. export imputed values, replace NAs in .xys files with imputed values (using knnreplace.kay_microarray.pl; or batch call KnnReplaceBatch2.pl)
#7. load modified .xys files into R using oligo, prune bad samples, run rma (this R file)
#8. use awk in shell to make datamatrix.txt (described in this R file)
#9. prune bad probes, rerun rma (this R file)?
#10. ready for analysis (this R file or as in Kay_analysis.txt or Kay_Final_analysis.txt - actually R code)

####Step 5####

#specify bioconductor source
source("http://bioconductor.org/biocLite.R")
#install packages
biocLite(c("affy", "oligo", "impute", "pdInfoBuilder"))
biocLite("limma")
biocLite("statmod")
biocLite("maanova")

#loaded libraries: oligo, limma
library(limma)
#library(affy)
library(oligo)
library(impute)
#library(genefilter)

#construct annotation library
library(pdInfoBuilder)
baseDir <- "~/GitHub/Cdiff_expression/remove bad features"
xys <- list.files(baseDir, pattern = ".xys", full.names = TRUE)[14] #choose xys file that doesn't need correction? Here 500952_Cycle10_532.xys
baseDir <- "~/GitHub/Cdiff_expression"
ndf <- list.files(baseDir, pattern = ".ndf", full.names = TRUE) #design file
seed <- new("NgsExpressionPDInfoPkgSeed", ndfFile = ndf, xysFile = xys, biocViews = "AnnotationData", organism = "knapweed", species = "Centaurea diffusa")
makePdInfoPackage(seed, destDir = baseDir)
#now install package. there is no need to load it in after. R will recognize that it's installed once you try to read the xys files
install.packages("~/GitHub/Cdiff_expression/pd.110405.cdiffusa.lz.exp/", repos=NULL, type="source")

#read in the data 
xys.files<-list.xysfiles("~/GitHub/Cdiff_expression/remove bad features", full.names=TRUE)
Cdifexpr<-read.xysfiles(xys.files)

#remove random probes	#identifies probes which only have NA values across every array
data<-exprs(Cdifexpr)
# data2<-data[apply(data,1,function(x)any(!is.na(x))),]
# data3<-data[apply(data,1,function(x)any(is.na(x))),]
data2<-data[rowSums(is.na(data)) < 108, ,drop=FALSE] #number here should match the number of samples (columns) in the data frame
#look at raw data
hist(data2, col = darkColors(132), lty = 1) #, xlim = c(100, 5000)


#k nearest neighbor imputation
knnout<-impute.knn(data2) #imputing with default values #k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=362436069
allkn<-knnout$data

#print out and replace original values with imputed values using kay's script (knnreplace.kay_mircroarray.pl)
write.table(allkn, file="~/GitHub/Cdiff_expression/remove bad features/knn_replaced.txt", quote=F, sep="\t")
#file too large to push to GitHub; stored locally in C:\Users\Kat\Documents\grad work\gene expression\remove bad features
#also, must add tab at beginning of file to align column headers properly

####Step 7####

#read in KNN imputed XYS files
xys.files<-list.xysfiles("~/GitHub/Cdiff_expression/knn_xys", full.names=TRUE)
Cdifxys<-read.xysfiles(xys.files)

#RMA - background correction, normalization, summarization of probes
ppData.Cdifxys=rma(Cdifxys)

#RMA is complete, now extract expression data
exprs.Cdifxys=exprs(ppData.Cdifxys)

#not sure these work properly
write.table(exprs.Cdifxys, file="Cdif_exprs_matrix.txt")
readtest<- read.delim("Cdif_exprs_matrix.txt", header=TRUE, sep = "\t")

####Step 8####
#extracting probe expression values from sunflower probe info files #DONE IN SHELL/CYGWIN NOT R
#awk command as below, written out with appropriate file names in probelist.txt
#awk '{a[FNR]=a[FNR] FS $8;t=(FNR>T)?FNR:t}END {for (i=1;i<=t;i++) print a[i]}' blah.probe blahblah.probe etc.probe > datamatrix.txt

####Step 9####
probes=read.delim("C:/Users/Kat/Documents/grad work/gene expression/34029_20120220_Knapweeds/probe_reports/datamatrix.txt", sep=" ", header=T, row.names=NULL)

#above output doesn't contain probe names, get them from a single probe file
a=read.delim("C:/Users/Kat/Documents/grad work/gene expression/34029_20120220_Knapweeds/probe_reports/509196_Cycle7_532.probe", sep="\t", header=T, row.names=NULL, skip=1) 
#duplicate probe ID names means we have to make them unique somehow
b=do.call(paste, c(a[c("PROBE_ID", "MATCH_INDEX")], sep = "_"))	
#concatenating the match index with the probe ID produces a human readable unique probe ID
b <- c(0,b)
row.names(probes)=b #assign unique row names to data matrix

#remove extra column and row introduced by awk. Remove value for 1st row, every column (which is "SEQ_URL")
probes=probes[-1,-1]

#REMOVE NON-RANDOM CONTROL PROBES FOR "BETTER" NORMALIZATION AND RMA
index1 <- with(probes, grepl("CPK", row.names(probes))) #which probes are CPK probes
index2 <- with(probes, grepl("AMPFS", row.names(probes))) #which probes are other non-random control probes
index3 <- with(probes, grepl("DELFS", row.names(probes))) #which probes are other non-random control probes
index4 <- with(probes, grepl("XENOTRACK", row.names(probes))) #which probes are other non-random control probes
index5 <- with(probes, grepl("EMPTY", row.names(probes))) #which probes are other non-random control probes
index6 <- with(probes, grepl("Syn[A-Z]_", row.names(probes))) #which probes are other non-random control probes
a=which(index1 == "FALSE") #I WANT TO REMOVE THESE
b=which(index2 == "FALSE") #I WANT TO REMOVE THESE
c=which(index3 == "FALSE") #I WANT TO REMOVE THESE
d=which(index4 == "FALSE") #I WANT TO REMOVE THESE
e=which(index5 == "FALSE") #I WANT TO REMOVE THESE
f=which(index6 == "FALSE") #I WANT TO REMOVE THESE
g=Reduce(intersect, list(a,b,c,d,e,f)) #only want those that are false in all cases
test=probes[g,]
test <- data.matrix(test, rownames.force=TRUE)

#RMA with limma this time. no using oligo without .xys inputs. 
#background correction
bgc=backgroundCorrect(test, method="normexp", verbose=T)
# bgc=backgroundCorrect(test, method="rma",verbose=T)
#quantile normalization
QNprobes=normalizeBetweenArrays(bgc, method="quantile")

#just gets the random probes on their own
#test=subset(QNsunprobes, grepl("RANDOM", row.names(QNsunprobes))) #just gets the random probes

#define random probes by index
QNprobes=as.data.frame(QNprobes)
row.names(QNprobes) <- row.names(test)
index1 <- with(QNprobes, grepl("RANDOM", row.names(QNprobes))) #which probes are RANDOM probes
a=which(index1 == "FALSE") #these are the non-random probes
b=which(index1 == "TRUE") #these are the random probes
rand=QNprobes[b,]
data=QNprobes[a,]

#define cutoff for 'unexpressed' <= 2 * (Standard deviation of random probes) + (mean random probes)
randmat <- data.matrix(rand, rownames.force=TRUE)
mean(randmat)
[1] 2876.698
sd(randmat)
[1] 1067.049
(1067.049 * 2) + 2876.698
[1] 5010.796	#remove genes with expression below this level

min(exprs.Cdifxys)
[1] 2.403344	
mean(exprs.Cdifxys)

#remove unexpressed


####Step 10######
#####using MAANOVA package to analyze data w a mixed model
#building this from MAANOVA tutorial

# Read hybridization design
design=read.table("C:/Users/Kat/Documents/GitHub/Cdiff_expression/experimentaldesign.txt", header=T, sep="\t") 
#expt design table includes expression and DNA arrays
"arrayID"  "SampleID" "Pop"      "Origin"   "Trt"      "Exp"      "TrtPool"  "Tmpt"
#arrayID needs to be Array
colnames(design)[1] <- "Array"

exprs.design <- droplevels(subset(design, Exp=="exprs"))

#load maanova library
library(maanova)

# Create a madata object #data is already normalized, background corrected, summarized. don't log transform.
madata=read.madata(exprs.Cdifxys, exprs.design, log.trans=F)

###############################################################################################################
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
