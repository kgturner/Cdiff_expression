#C.diffusa expression - clustering
#2/16/15

####kmeans clustering####
#pca of all genes
#t0
library("ggplot2")
library("grid") 
library("gridBase")

###load data, if neccessary####
#read data table, probably slow
#on cluster
exprs.df <- read.table("~/Centaurea_diffusa_expression/Cdifexprs_lme4dat.txt", header=T, sep="\t") 
#local
exprs.df <- read.table("C:/Users/Kat/Documents/GitHub/Cdiff_expression/Cdifexprs_lme4dat.txt", header=T, sep="\t") 
#test!
test <- read.table("test_lme4dat.txt", header=T, sep="\t") 
testc <- subset(test, Trt=="control")
testc <- droplevels(testc)

testd <- subset(test, Trt=="drought")

#control test
test.pca <- prcomp(test[c(15:20)], center=TRUE, retx=T, scale.=TRUE)
summary(test.pca)
# Importance of components:
#   PC1    PC2    PC3    PC4     PC5     PC6
# Standard deviation     1.3479 1.2079 1.1107 0.8888 0.63962 0.53972
# Proportion of Variance 0.3028 0.2432 0.2056 0.1317 0.06819 0.04855
# Cumulative Proportion  0.3028 0.5460 0.7516 0.8833 0.95145 1.00000
 
#visualize components
plot(test.pca, main="(a) Screeplot, All Occurrences", xlab="Principal component", ylim=c(0,7))

# biplot(test.pca)
# biplot(test.pca,  main="PCA analysis of gene expr data", choices=c(1,3))

#see bottom for figure
# 
# variances of the principal components:
apply(test.pca$x, 2, var)
# PC1       PC2       PC3       PC4       PC5       PC6 
# 1.8168919 1.4590836 1.2335725 0.7900455 0.4091118 0.2912946 

loadings <- test.pca$rotation[,1]
sort(abs(loadings), decreasing=TRUE)
# Contig1 Contig10000  Contig1000   Contig100 Contig10001    Contig10 
# 0.6469863   0.5014620   0.4492758   0.3193297   0.1361187   0.0871555 

#find top loadings (for PC2)
loadings2 <- test.pca$rotation[,2]
sort(abs(loadings2), decreasing=TRUE)
# Contig10001        Tmpt  Contig1000         PC1   Contig100    Contig10     Contig1 
# 0.62045772  0.50808524  0.46746503  0.28090670  0.14978544  0.13680545  0.11931250 
# Contig10000 
# 0.06373843 

# #find top loadings (for PC3)
# loadings3 <- allclim.pca$rotation[,3]
# sort(abs(loadings3), decreasing=TRUE)

#proportional contributions of each bioclim to each PC
#If you want this as a relative contribution then sum up the loadings per column and 
#express each loading as a proportion of the column (loading) sum, taking care to use 
#the absolute values to account for negative loadings.

sweep(abs(test.pca$rotation),2, colSums(abs(test.pca$rotation)),"/")
#             PC1        PC2        PC3        PC4        PC5        PC6        PC7
# Tmpt        0.21652376 0.08968505 0.10807235 0.11515466 0.06978513 0.17165913 0.24804065
# PC1         0.11971018 0.05378969 0.18679140 0.04867560 0.34922940 0.01062098 0.08524038
# Contig1     0.05084578 0.24238612 0.07059099 0.02118569 0.06910217 0.27095543 0.17181873
# Contig10    0.05830051 0.10473799 0.22473090 0.10448648 0.21718817 0.17968800 0.05171808
# Contig100   0.06383202 0.13351128 0.11057210 0.39524371 0.06792274 0.03885579 0.10167091
# Contig1000  0.19921320 0.15116696 0.09920496 0.07067838 0.09187255 0.01887156 0.12681840
# Contig10000 0.02716254 0.18794147 0.15022528 0.20632834 0.09467758 0.27957711 0.02026756
# Contig10001 0.26441201 0.03678144 0.04981202 0.03824714 0.04022226 0.02977200 0.19442529

# # #get top 4 PCs
# PC1 <- as.matrix(allclim.pca$x[,1])
# PC2 <- as.matrix(allclim.pca$x[,2])
# PC3 <- as.matrix(allclim.pca$x[,3])
# # # PC4 <- as.matrix(Frclim.pca$x[,4])
# 
# allclim2 <- cbind(allclim, PC1, PC2, PC3)
# # 
# # #write table
# write.table(allclim2, file="Cdif_allocc_bioclimPCA.txt")
# # # 
# allclim2 <- read.table("Cdif_allocc_bioclimPCA.txt", header=TRUE)
# 
# ####all occ main fig; 95% conf limits of clusters####
# # http://stackoverflow.com/questions/20260434/test-significance-of-clusters-on-a-pca-plot
# # draw 95% confidence ellipses around clusters. Note that stat_ellipse(...) uses the bivariate t-distribution.
scores <- test.pca$x[,1:3]                        # scores for first three PC's

# k-means clustering [assume 6 clusters] for all data in one plot
km     <- kmeans(scores, centers=6, nstart=10)
ggdata <- data.frame(scores, Cluster=km$cluster, Origin=test$Origin, Trt=test$Trt,
                     Pop=test$Pop, PopTrtPool=test$PopTrtPool, Tmpt=test$Tmpt)
levels(ggdata$Origin)[levels(ggdata$Origin)=="inv"] <- "Invasive C. diffusa"
levels(ggdata$Origin)[levels(ggdata$Origin)=="nat"] <- "Native C. diffusa"

# stat_ellipse is not part of the base ggplot package
source("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R") 

#centroid based on origin
centroids <- aggregate(cbind(PC1,PC2)~Cluster+Origin+Pop,data=ggdata,mean)

#PC1 vs PC2, blind cluster
plot <- ggplot(ggdata, aes_string(x="PC1", y="PC2")) +
  geom_point(aes(color=factor(Pop),shape=Origin), size=3) +
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
               geom="polygon", level=0.95, alpha=0.2) #+
#   guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) #+
#   facet_grid(Trt ~ Tmpt)
plot

# k-means clustering [assume 36 clusters] for faceted data
km     <- kmeans(scores, centers=36, nstart=40)
ggdata <- data.frame(scores, Cluster=km$cluster, Origin=test$Origin, Trt=test$Trt,
                     Pop=test$Pop, PopTrtPool=test$PopTrtPool, Tmpt=test$Tmpt)
levels(ggdata$Origin)[levels(ggdata$Origin)=="inv"] <- "Invasive C. diffusa"
levels(ggdata$Origin)[levels(ggdata$Origin)=="nat"] <- "Native C. diffusa"

# stat_ellipse is not part of the base ggplot package
source("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R") 

#centroid based on origin
centroids <- aggregate(cbind(PC1,PC2)~Origin+Pop+Trt+Tmpt,data=ggdata,mean)

#95% plot, by pop
Oplot <- ggplot(ggdata, aes_string(x="PC1", y="PC2")) +
  geom_point(aes(color=factor(Pop),shape=Origin), size=3) +
  guides(color=guide_legend("Pop"),fill=guide_legend("Pop"))+
#   stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
#                geom="polygon", level=0.95, alpha=0.2) +
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=Pop, shape=Origin), size=8) +
  #coord_cartesian(ylim = c(-6.5, 8.5)) + 
  facet_grid(Trt ~ Tmpt) +
  theme_bw() + 
  theme(legend.justification=c(1,0), legend.position=c(1,0),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size = 10))

Oplot
# ggsave("KTurnerFig4.pdf", width=6.65, height = 5)
# ggsave("KTurnerFig4.png", width=6.65, height = 5)
# 
# svg("KTurnerFig4.svg", width=6.65, height=5, pointsize = 12)
# Oplot
# dev.off()





####trial####

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
