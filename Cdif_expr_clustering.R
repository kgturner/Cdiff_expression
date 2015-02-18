#C.diffusa expression - clustering
#2/16/15
#based on Gillespie et al., 2010. Github repo at https://github.com/csgillespie/bmc-microarray
#and http://stats.stackexchange.com/questions/3271/clustering-genes-in-a-time-course-experiment

source("http://bioconductor.org/biocLite.R")
biocLite(c("hopach", "pvclust","Mfuzz", "GeneNet"))

library("hopach")
# library("Biobase")
library("pvclust")
library("gplots")
library("Mfuzz")
library("GeneNet")

#load data
PC1q_intsigdf <- read.table("PC1_sigint_df.txt", header=T)
# PC1q_Osigdf <- read.table("PC1_sigOrigin_df.txt", header=T)
# PC1q_sigdf <- read.table("PC1_sigPC1_df.txt", header=T)
# PC1q_trtsigdf <- read.table("PC1_sigtrt_df.txt", header=T)
# #pop means?
# popintdf <- read.table("PC1_sigint_popMeans.txt", header=T)

# #for Mfuzz, subset by trt and origin and average over pools and pops
# # #example
# library(plyr)
# test <- PC1q_intsigdf[,c(1:241)]
# test2 <- subset(test, Origin=="nat"&Trt=="control")
# test3 <- ddply(test2, .(Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))
# test3[,1:30]
# 
# test3M <- as.matrix(t(test3[,c(9:235)]))
# all(colnames(t(test3[,c(9:235)]))==row.names(test3))
# test3eset <- new("ExpressionSet",exprs = test3M) #phenoData = as(test3[,1:8], "AnnotatedDataFrame")
# test3eset_st <- standardise(test3eset)
# cl = mfuzz(test3eset_st, c = 8, m = 1.25) #c=# of clusters
# mfuzz.plot(test3eset_st, cl = cl, mfrow = c(2, 4), new.window = FALSE)

####pick number of clusters?####
library(hopach)
#example for nat control, sig int
#distance matrix
gene.dist<-distancematrix(test3eset_st,"cosangle")
dim(gene.dist)
gene.hobj<-hopach(test3eset_st,dmat=gene.dist)
gene.hobj$clust$k
[1] 32
table(gene.hobj$clust$labels)
gene.hobj$clust$sizes
dplot(gene.dist,gene.hobj,ord="final",main="Gene Distance Matrix",showclusters=FALSE)

#nat dr, sig int
gene.dist<-distancematrix(natdr3eset_st,"cosangle")
dim(gene.dist)
gene.hobj<-hopach(natdr3eset_st,dmat=gene.dist)
gene.hobj$clust$k
[1] 43
table(gene.hobj$clust$labels)
gene.hobj$clust$sizes
dplot(gene.dist,gene.hobj,ord="final",main="Gene Distance Matrix",showclusters=FALSE)

#inv ctrl, sig int
gene.dist<-distancematrix(invctrl3eset_st,"cosangle")
dim(gene.dist)
gene.hobj<-hopach(invctrl3eset_st,dmat=gene.dist)
gene.hobj$clust$k
[1] 20
table(gene.hobj$clust$labels)
gene.hobj$clust$sizes
dplot(gene.dist,gene.hobj,ord="final",main="Gene Distance Matrix",showclusters=FALSE)

#inv dr, sig int
gene.dist<-distancematrix(invdr3eset_st,"cosangle")
dim(gene.dist)
gene.hobj<-hopach(invdr3eset_st,dmat=gene.dist)
gene.hobj$clust$k
[1]89
table(gene.hobj$clust$labels)
gene.hobj$clust$sizes
dplot(gene.dist,gene.hobj,ord="final",main="Gene Distance Matrix",showclusters=FALSE)


####soft clustering with Mfuzz####
#necessary to standardize???
# standardise their measurements by taking the expression level of the mutant strain 
# (at each timepoint) relative to the wild-type at time t = 0
c_probe_data = yeast.matrix[ii, ]
# Average of WT
wt_means = apply(c_probe_data[, 16:30], 1, mean)
m = matrix(nrow = dim(c_probe_data)[1], ncol = 5)
for (i in 1:5) {
  mut_rep = c(i, i + 5, i + 10)
  m[, i] = rowMeans(c_probe_data[, mut_rep]) - wt_means
}
colnames(m) = sort(unique(exp_fac$tps))

#native control, sig Origin*Trt
library(plyr)
test <- PC1q_intsigdf[,c(1:241)]
test2 <- subset(test, Origin=="nat"&Trt=="control")
test3 <- ddply(test2, .(Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))
test3[,1:30]

test3M <- as.matrix(t(test3[,c(9:235)]))
all(colnames(t(test3[,c(9:235)]))==row.names(test3))
test3eset <- new("ExpressionSet",exprs = test3M) #phenoData = as(test3[,1:8], "AnnotatedDataFrame")
test3eset_st <- standardise(test3eset)
cl = mfuzz(test3eset_st, c = 2, m = 1.25) #c=# of clusters
mfuzz.plot(test3eset_st, cl = cl, mfrow = c(1, 2), new.window = FALSE)

#probes present within each cluster
cluster = 1
cl[[4]][, cluster]

#native, drought, sig int Origin*Trt
library(plyr)
natdr <- PC1q_intsigdf[,c(1:241)]
natdr2 <- subset(natdr, Origin=="nat"&Trt=="drought")
natdr3 <- ddply(natdr2, .(Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))
natdr3[,1:30]

natdr3M <- as.matrix(t(natdr3[,c(9:235)]))
all(colnames(t(natdr3[,c(9:235)]))==row.names(natdr3))
natdr3eset <- new("ExpressionSet",exprs = natdr3M) #phenoData = as(test3[,1:8], "AnnotatedDataFrame")
natdr3eset_st <- standardise(natdr3eset)
cl = mfuzz(natdr3eset_st, c = 12, m = 1.25) #c=# of clusters
mfuzz.plot(natdr3eset_st, cl = cl, mfrow = c(3, 4), new.window = FALSE)

#probes present within each cluster
cluster = 1
cl[[4]][, cluster]

#invasive, control, sig int Origin*Trt
library(plyr)
invctrl <- PC1q_intsigdf[,c(1:241)]
invctrl2 <- subset(invctrl, Origin=="inv"&Trt=="control")
invctrl3 <- ddply(invctrl2, .(Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))
invctrl3[,1:30]

invctrl3M <- as.matrix(t(invctrl3[,c(9:235)]))
all(colnames(t(invctrl3[,c(9:235)]))==row.names(invctrl3))
invctrl3eset <- new("ExpressionSet",exprs = invctrl3M) #phenoData = as(test3[,1:8], "AnnotatedDataFrame")
invctrl3eset_st <- standardise(invctrl3eset)
cl = mfuzz(invctrl3eset_st, c = 12, m = 1.25) #c=# of clusters
mfuzz.plot(invctrl3eset_st, cl = cl, mfrow = c(3, 4), new.window = FALSE)

#probes present within each cluster
cluster = 1
cl[[4]][, cluster]

#invasive, control, sig int Origin*Trt
library(plyr)
invdr <- PC1q_intsigdf[,c(1:241)]
invdr2 <- subset(invdr, Origin=="inv"&Trt=="drought")
invdr3 <- ddply(invdr2, .(Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))
invdr3[,1:30]

invdr3M <- as.matrix(t(invdr3[,c(9:235)]))
all(colnames(t(invdr3[,c(9:235)]))==row.names(invdr3))
invdr3eset <- new("ExpressionSet",exprs = invdr3M) #phenoData = as(test3[,1:8], "AnnotatedDataFrame")
invdr3eset_st <- standardise(invdr3eset)
cl = mfuzz(invdr3eset_st, c = 12, m = 1.25) #c=# of clusters
mfuzz.plot(invdr3eset_st, cl = cl, mfrow = c(3, 4), new.window = FALSE)

#probes present within each cluster
cluster = 1
cl[[4]][, cluster]


####regulatory network inference with GeneNet####
# data are stored in a matrix m, where the rows are genes and the fifthteen columns are the arrays. 
# we have time course data, so rearrange the row order according
# to the time points. The first three rows of the resulting matrix mnew are data on the three mutant
# arrays at time t = 0, the next three arrays at time point t = 60, and so on
# exp_fac = with(exp_fac, exp_fac[order(strain, tps, replicates), ])

genet1 <- with(PC1q_intsigdf, PC1q_intsigdf[order(Origin, Tmpt, TrtPool),])
genet1 <- subset(genet1, Origin=="nat"&Trt=="control")

# Construct a longitudinal object
library(GeneNet)
genetM <- as.matrix(t(genet1[,c(15:241)]))
mlong = as.longitudinal(t(genetM), repeats = 9, time = 0:2)
# Compute partial correlations
pcor.dyn = ggm.estimate.pcor(mlong, method = "dynamic")
# Assign (local) fdr values to all possible edges
m.edges = network.test.edges(pcor.dyn, direct = TRUE)
# Construct graph containing top edges
m.net = extract.network(m.edges, method.ggm = "number", cutoff.ggm = 100)
# Construct a Graphviz dot file
rnames = vector("list", length(1))
rnames = rownames(genetM)
network.make.dot(filename = "net.dot", m.net, rnames, main = "nat, control, sig int, Network")

igr = network.make.igraph(m.net, rnames)
plot(igr, main="nat, control, sig int, Network", layout=layout.fruchterman.reingold,
     edge.arrow.size=0.5, vertex.size=9, vertex.label.cex=0.7)

#replace contig names with TAIR accessions?
