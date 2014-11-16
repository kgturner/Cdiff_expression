#Expression/Maternal effects comparison_PCA and climate data


#REML, using lme4
#mixed effect models 
# library(plyr)
# library(lme4.0)
# library(lsmeans)
# library("ggplot2")
# library("grid") 
# library("gridBase")

####PCA fun times####
Mfclim <- read.table("MfbioclimPCAdat.txt", header=TRUE)
exprs.clim <- subset(Mfclim, Pop%in%c("BG001", "TR001","RU008","CA001", "US001", "US003"), select=1:24)

#Origin and longitude artificially separates groups...
# Mfclim.pca <- prcomp(Mfbioclim[c(2:22)], center=TRUE, scale=TRUE)
exprs.clim.pca <- prcomp(exprs.clim[c(2:22)], center=TRUE, scale=TRUE)
summary(exprs.clim.pca)
# Importance of components:
#                           PC1    PC2    PC3     PC4     PC5       PC6
# Standard deviation     3.0179 2.3410 2.0838 1.15953 0.85161 7.932e-16
# Proportion of Variance 0.4337 0.2610 0.2068 0.06402 0.03453 0.000e+00
# Cumulative Proportion  0.4337 0.6947 0.9014 0.96547 1.00000 1.000e+00

#visualize components
plot(exprs.clim.pca, main="Variances of each principle component of climate", xlab="Principal component", ylim=c(0,7))
# screeplot(Frclim.pca, type="lines")
biplot(exprs.clim.pca)

#see bottom for figure

# variances of the principal components:
apply(exprs.clim.pca$x, 2, var)
#         PC1          PC2          PC3          PC4          PC5          PC6 
# 9.107895e+00 5.480248e+00 4.342122e+00 1.344500e+00 7.252347e-01 3.526378e-32 

# # biplot(Mfclim.pca, var.axes=FALSE, main="PCA analysis of climate data")
# biplot(exprs.clim.pca, var.axes=TRUE, main="PCA analysis of climate data", cex=c(1,2), col=c(Frclimdat$Origin,"red"))
# subset(Mfclim, Pop%in%c("SERG","CA001","US001","US002"), select=c(bio5,bio17,bio14,bio12))
# subset(Mfclim, Pop%in%c("RU008","SERG","TR001","CA001"), select=c(bio19,bio6,bio4,bio11))
# biplot(Mfclim.pca, var.axes=FALSE, main="PCA analysis of climate data", choices=c(1,3))
# subset(Mfclim, Pop%in%c("US001","US003","GR002","BG001"), select=c(bio1,bio3,bio15,bio2))
#see bottom for figure

#get top 3 PCs
PC1 <- as.matrix(exprs.clim.pca$x[,1])
PC2 <- as.matrix(exprs.clim.pca$x[,2])
PC3 <- as.matrix(exprs.clim.pca$x[,3])

exprs.clim <- cbind(exprs.clim, PC1, PC2, PC3)

#find top loadings (for PC1)
loadings <- exprs.clim.pca$rotation[,1]
sort(abs(loadings), decreasing=TRUE)
#       bio7       bio19        bio4       bio12       bio16       bio13        bio5        bio6         alt       bio17       bio11       bio14 
# 0.316560371 0.307971387 0.302076990 0.294497438 0.282495007 0.279244047 0.268251752 0.258418562 0.231386476 0.228779880 0.227676373 0.203608993 
# bio10        bio2    Latitude        bio3        bio9        bio8        bio1       bio18       bio15 
# 0.198324789 0.193318619 0.185396832 0.087759295 0.070945421 0.066687446 0.050592759 0.021235734 0.001206422 
BIO7 = Temperature Annual Range (BIO5-BIO6)
BIO19 = Precipitation of Coldest Quarter
BIO4 = Temperature Seasonality (standard deviation *100)
BIO12 = Annual Precipitation

#find top loadings (for PC2)
loadings2 <- exprs.clim.pca$rotation[,2]
sort(abs(loadings2), decreasing=TRUE)
#       bio3       bio2       bio8       bio9      bio18      bio15      bio14      bio17   Latitude      bio10       bio1       bio4      bio19 
# 0.40599722 0.33517677 0.31644489 0.29970539 0.29609469 0.28933246 0.28753175 0.25225216 0.24976603 0.22013759 0.15477109 0.14775442 0.14013676 
# bio5        alt       bio7      bio16      bio13      bio12       bio6      bio11 
# 0.12295028 0.11544500 0.08084751 0.05151160 0.02345269 0.02247143 0.01325613 0.00623678 
BIO3 = Isothermality (BIO2/BIO7) (* 100)
BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
BIO8 = Mean Temperature of Wettest Quarter
BIO9 = Mean Temperature of Driest Quarter

#find top loadings (for PC3)
loadings3 <- exprs.clim.pca$rotation[,3]
sort(abs(loadings3), decreasing=TRUE)
# bio1       bio11       bio18        bio6        bio9       bio10        bio8    Latitude       bio13       bio16        bio5         alt 
# 0.424181557 0.342687982 0.321366393 0.295437537 0.288419778 0.270379804 0.234966770 0.233497494 0.217498824 0.207345523 0.207217645 0.202523882 
# bio15       bio12        bio4        bio7       bio14        bio2       bio19       bio17        bio3 
# 0.187274951 0.155343910 0.089744774 0.057829287 0.028535784 0.016194890 0.010702512 0.002012817 0.001859683 
BIO1 = Annual Mean Temperature
BIO11 = Mean Temperature of Coldest Quarter
BIO18 = Precipitation of Warmest Quarter
BIO6 = Min Temperature of Coldest Month

#proportional contributions of each bioclim to each PC
#If you want this as a relative contribution then sum up the loadings per column and 
#express each loading as a proportion of the column (loading) sum, taking care to use 
#the absolute values to account for negative loadings.

sweep(abs(exprs.clim.pca$rotation),2, colSums(abs(exprs.clim.pca$rotation)),"/")
#                   PC1         PC2          PC3         PC4         PC5         PC6
# alt      0.0567619745 0.030132286 0.0533656703 0.120768559 0.042267011 0.011099892
# bio1     0.0124110317 0.040396783 0.1117731544 0.040702304 0.055555802 0.060592568
# bio10    0.0486515323 0.057458086 0.0712459160 0.020913573 0.065682263 0.017519819
# bio11    0.0558518400 0.001627861 0.0902993447 0.031442196 0.011158624 0.057817586
# bio12    0.0722438767 0.005865265 0.0409336015 0.046203116 0.075109700 0.070963298
# bio13    0.0685020305 0.006121382 0.0573116139 0.016608243 0.080179766 0.022138760
# bio14    0.0499478130 0.075048627 0.0075192676 0.099176707 0.018366936 0.021912223
# bio15    0.0002959502 0.075518631 0.0493475298 0.151062654 0.033747208 0.029899147
# bio16    0.0692995315 0.013445036 0.0546361877 0.008723885 0.076906622 0.031510303
# bio17    0.0561225441 0.065840306 0.0005303835 0.102996802 0.002636005 0.016136493
# bio18    0.0052093892 0.077283639 0.0846810400 0.008158576 0.074253242 0.025692882
# bio19    0.0755492036 0.036577079 0.0028201450 0.030298693 0.032892503 0.021455329
# bio2     0.0474234566 0.087484448 0.0042674037 0.040855039 0.035847063 0.091488399
# bio3     0.0215284444 0.105969283 0.0004900322 0.007612543 0.046241813 0.012776914
# bio4     0.0741032350 0.038565364 0.0236480260 0.007014427 0.033580976 0.150379295
# bio5     0.0658054841 0.032091237 0.0546024913 0.041034734 0.063122572 0.180176236
# bio6     0.0633932806 0.003459981 0.0778487063 0.025074749 0.011107484 0.009308032
# bio7     0.0776561879 0.021101998 0.0152381962 0.039595482 0.030566335 0.037582565
# bio8     0.0163592581 0.082595241 0.0619144718 0.091829261 0.051986945 0.025511569
# bio9     0.0174037922 0.078226066 0.0759995050 0.055087041 0.065189821 0.052364206
# Latitude 0.0454801440 0.065191401 0.0615273131 0.014841416 0.093601309 0.053674484

#write table
write.table(exprs.clim, file="CdifExprs.BioclimPCAdat.txt")
# 
exprs.clim <- read.table("CdifExprs.BioclimPCAdat.txt", header=TRUE)

# ##########figures#######
# #nat and inv and sk colors:"#F8766D","#00BFC4"
# library("ggplot2")
# library("grid") 
# library(gridBase)
# 
# ###PC1 vs PC2####
# #pts instead of labels for pops
# data <- data.frame(obsnames=row.names(Mfclim.pca$x), Mfclim.pca$x)
# data <- merge(data, Mfclim[c(1,24)],by.x="obsnames", by.y="Pop")
# levels(data$Origin)[levels(data$Origin)=="inv"] <- "Invasive C. diffusa"
# levels(data$Origin)[levels(data$Origin)=="nat"] <- "Native C. diffusa"
# levels(data$Origin)[levels(data$Origin)=="sk"] <- "Native C. stoebe"
# 
# # pdf("KTurnerFig2.pdf", useDingbats=FALSE, width=13.38)
# png("MfClimatePC1vPC2.png",width=800, height = 600, pointsize = 16)
# # postscript("KTurnerFig2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 13.38)
# 
# pPC2 <- ggplot(data, aes_string(x="PC1", y="PC2")) + 
#   geom_point(aes(shape=Origin, color=Origin), size=3) +
#   #   scale_x_continuous(expand = c(0,1)) #+
#   theme_bw()+
#   theme(legend.justification=c(1,0), legend.position=c(1,0))
# 
# # plot
# 
# pPC2 <- pPC2 + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
# datapc <- data.frame(varnames=rownames(Mfclim.pca$rotation), Mfclim.pca$rotation)
# mult <- min(
#   (max(data[,"PC2"]) - min(data[,"PC2"])/(max(datapc[,"PC2"])-min(datapc[,"PC2"]))),
#   (max(data[,"PC1"]) - min(data[,"PC1"])/(max(datapc[,"PC1"])-min(datapc[,"PC1"])))
# )
# datapc <- transform(datapc,
#                     v1 = .7 * mult * (get("PC1")),
#                     v2 = .7 * mult * (get("PC2"))
# )
# 
# pPC2 <- pPC2 + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), 
#                                          size = 6, vjust=1, color="gray47", alpha=0.75)
# pPC2 <- pPC2 + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), 
#                             arrow=arrow(length=unit(0.2,"cm")), alpha=0.4, color="gray47")+
#   ggtitle("(b)")+theme(plot.title = element_text(lineheight=2, face="bold"))
# pPC2
# dev.off()
# 
# ####PC1 vs PC3####
# png("MfClimatePCA1v3.png",width=800, height = 600, pointsize = 16)
# # postscript("KTurnerFig2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 13.38)
# 
# pPC3 <- ggplot(data, aes_string(x="PC1", y="PC3")) + 
#   geom_point(aes(shape=Origin, color=Origin), size=3) +
#   #   scale_x_continuous(expand = c(0,1)) #+
#   theme_bw()+
#   theme(legend.position="none")
# #   theme(legend.justification=c(1,0), legend.position=c(1,0))
# 
# # plot
# 
# pPC3 <- pPC3 + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
# 
# datapc <- data.frame(varnames=rownames(Mfclim.pca$rotation), Mfclim.pca$rotation)
# mult <- min(
#   (max(data[,"PC3"]) - min(data[,"PC3"])/(max(datapc[,"PC3"])-min(datapc[,"PC3"]))),
#   (max(data[,"PC1"]) - min(data[,"PC1"])/(max(datapc[,"PC1"])-min(datapc[,"PC1"])))
# )
# datapc <- transform(datapc,
#                     v1 = .7 * mult * (get("PC1")),
#                     v2 = .7 * mult * (get("PC3"))
# )
# 
# pPC3 <- pPC3  +coord_equal(ratio=6.1/4)+ geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), 
#                                                    size = 6, vjust=1, color="gray47", alpha=0.75)
# 
# pPC3 <- pPC3 + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), 
#                             arrow=arrow(length=unit(0.2,"cm")), alpha=0.4, color="gray47")+
#   ggtitle("(c)")+theme(plot.title = element_text(lineheight=2, face="bold"))
# pPC3
# dev.off()
# 
# 
# ####supp. fig multiplot####
# # pdf("KTurnerFig2.pdf", useDingbats=FALSE, width=13.38)
# png("KTurnerSup_VanPCA.png",width=800, height = 800, pointsize = 12)
# # postscript("KTurnerFig2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 13.38)
# 
# gl <- grid.layout(2, 2)
# #                   widths=unit(c(4, 4), "inches"),
# #                   heights=unit(c(4, 4), "inches"))
# pushViewport(viewport(layout=gl))
# 
# par(mfcol=c(2,2))
# pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
# plot(Mfclim.pca, main="(a)", 
#      xlab="Principal component", ylim=c(0,8), cex.main=1.3)
# popViewport()
# 
# plot.new()
# pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
# vp1 <- plotViewport(c(0,0,0,0))
# print(pPC2, vp=vp1)
# popViewport()
# 
# plot.new()
# pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
# vp1 <- plotViewport(c(0,0,0,0))
# print(pPC3, vp=vp1)
# 
# 
# dev.off()