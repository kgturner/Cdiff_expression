#C.diffusa expression - differential expression
#1/19/15

####look at results####
#local
PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") 
# PC1qdr <- read.table("lme4_qval_PC1_dr.txt", header=T, sep="\t") 
# Latq <- read.table("lme4_qval_lat.txt", header=T, sep="\t")
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
PC1q_trtsig <- subset(PC1q, trtQsig==TRUE&intQsig==FALSE)#contigs with sig trt
# PC1q_trtsig <- subset(PC1q, trtQsig==TRUE)

pc1intV <- as.vector(PC1q_intsig$Contig)
PC1q_intsigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1intV)
PC1q_intsigdf <-cbind(exprs.df[,1:14], PC1q_intsigdf)

pc1oV <- as.vector(PC1q_Osig$Contig)
PC1q_Osigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1oV)
PC1q_Osigdf <-cbind(exprs.df[,1:14], PC1q_Osigdf)

pc1trtV <- as.vector(PC1q_trtsig$Contig)
#
PC1q_trtsigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1trtV)
PC1q_trtsigdf <-cbind(exprs.df[,1:14], PC1q_trtsigdf)

pc1V <- as.vector(PC1q_sig$Contig)
PC1q_sigdf <- subset(exprs.df, select=colnames(exprs.df)%in%pc1V)
PC1q_sigdf <-cbind(exprs.df[,1:14], PC1q_sigdf)

# #do I need set where origin AND trt sig?
# intersect(PC1q_Osigdf, PC1q_trtsigdf

write.table(PC1q_intsigdf, file="PC1_sigint_df.txt", sep="\t")
write.table(PC1q_Osigdf, file="PC1_sigOrigin_df.txt", sep="\t")
write.table(PC1q_sigdf, file="PC1_sigPC1_df.txt", sep="\t")
write.table(PC1q_trtsigdf, file="PC1_sigtrt_df.txt", sep="\t")

# #lat, not run
# Latq_intsig <- subset(Latq,intQsig==TRUE )
# Latq_Osig <- subset(Latq,originQsig==TRUE&intQsig==FALSE )
# Latq_sig <- subset(Latq, covQsig==TRUE)
#not working
# Latq_intsigdf <- subset(exprs.df, Contig%in%Latq_intsig$Contig)
# Latq_Osigdf <- subset(exprs.df, Contig%in%Latq_Osig$Contig)
# Latq_sigdf <- subset(exprs.df, Contig%in%Latq_sig$Contig)
# write.table(Latq_intsigdf, file="Lat_sigint_df.txt", sep="\t")
# write.table(Latq_Osigdf, file="Lat_sigOrigin_df.txt", sep="\t")
# write.table(Latq_sigdf, file="Lat_sigPC1_df.txt", sep="\t")


summary(PC1q_Osig)
# summary(Latq_Osig)


####pop means####
library(plyr)
# test <- PC1q_intsigdf[,c(1:34)]
# test2 <- ddply(test, .(Pop, Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))

popintdf <- ddply(PC1q_intsigdf, .(Pop, Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))
popintdf$OriginTrt <- as.factor(paste0(popintdf$Origin, "_", popintdf$Trt))
popintdf$PopTrtTmpt <- as.factor(paste0(popintdf$Pop, "_", popintdf$Trt,"_",popintdf$Tmpt))

write.table(popintdf, file="PC1_sigint_popMeans.txt", sep="\t")
popintdf <- read.table("PC1_sigint_popMeans.txt", header=T)

popOdf <- ddply(PC1q_Osigdf, .(Pop, Origin, Trt, Tmpt),numcolwise(mean, na.rm=TRUE))
popOdf$OriginTrt <- as.factor(paste0(popOdf$Origin, "_", popOdf$Trt))
popOdf$PopTrtTmpt <- as.factor(paste0(popOdf$Pop, "_", popOdf$Trt,"_",popOdf$Tmpt))

write.table(popOdf, file="PC1_sigOrigin_popMeans.txt", sep="\t")
popOdf <- read.table("PC1_sigOrigin_popMeans.txt", header=T)


####up and down expr between ranges tmpt 2, sig int####
####int using pop means, tmpt 2,control####
popintdf2c <- subset(popintdf, Tmpt==2&Trt=="control", select=c(1,2,10:236))
head(popintdf2c)

library(reshape2)
testdir <- reshape(popintdf2c,direction="long", varying=list(ExprVal=c(3:229)), times=colnames(popintdf2c[3:229]))
head(testdir)

library(plyr)
testdir5 <- ddply(testdir, .(time,Origin), summarise, mean(Contig1007))
head(testdir5)

testdir6 <- reshape(testdir5, direction="wide", idvar="time", timevar="Origin")
head(testdir6)
colnames(testdir6)[2] <- "InvExprVal"
colnames(testdir6)[3] <- "NatExprVal"
colnames(testdir6)[1] <- "Contig"
testdir6$InvUp <- testdir6$InvExprVal > testdir6$NatExprVal

int2csummary <- testdir6
summary(int2csummary)
write.table(int2csummary, file="sigint_popMeans_sumT2Control.txt", sep="\t")

####int using pop means, tmpt 2,drought####
popintdf2d <- subset(popintdf, Tmpt==2&Trt=="drought", select=c(1,2,10:236)) #last 2 cols are compound names, i.e. OriginTrt
head(popintdf2d[,1:6])
testdir <- reshape(popintdf2d,direction="long", varying=list(3:229), times=colnames(popintdf2c[3:229]))
head(testdir)
testdir5 <- ddply(testdir, .(time,Origin), summarise, mean(Contig1007))
head(testdir5)
testdir6 <- reshape(testdir5, direction="wide", idvar="time", timevar="Origin")
head(testdir6)
colnames(testdir6)[2] <- "InvExprVal"
colnames(testdir6)[3] <- "NatExprVal"
colnames(testdir6)[1] <- "Contig"
testdir6$InvUp <- testdir6$InvExprVal > testdir6$NatExprVal

int2dsummary <- testdir6
head(int2dsummary)
summary(int2dsummary)
write.table(int2dsummary, file="sigint_popMeans_sumT2drought.txt", sep="\t")

####allowing for nearly equal things, sig int####
#floating points in comparisons... don't do equals, just think of less than
close_enough <- function(x, y, tolerance=sqrt(.Machine$double.eps)) {
  abs(x - y) <= tolerance
}

int2dsummary <- read.table(file="sigint_popMeans_sumT2drought.txt", header=T)
int2dsummaryNE <- int2dsummary
int2dsummaryNE$OriginsEq <- close_enough(int2dsummaryNE$InvExprVal, int2dsummaryNE$NatExprVal, tolerance=0.1)
int2dsummaryNE$InvDefUp <- int2dsummaryNE$InvUp 
int2dsummaryNE[int2dsummaryNE$OriginsEq==TRUE,]$InvDefUp <- FALSE
int2dsummaryNE$InvDefDn <- int2dsummaryNE$InvUp#FALSE means down
int2dsummaryNE[int2dsummaryNE$OriginsEq==TRUE,]$InvDefDn <- TRUE
write.table(int2dsummaryNE, file="sigint_popMeans_sumT2droughtNE.txt", sep="\t")

#####sig Origin using pop means, time point 0, both trt (equivalent)####
popOdf0 <- subset(popOdf, Tmpt==0, select=c(1,2,10:594))
popOdf0[1:6,1:20]
testdir <- reshape(popOdf0,direction="long", varying=list(ExprVal=c(3:587)), times=colnames(popOdf0[3:587]))
head(testdir)
testdir5 <- ddply(testdir, .(time,Origin), summarise, mean(Contig10161))
head(testdir5)
testdir6 <- reshape(testdir5, direction="wide", idvar="time", timevar="Origin")
head(testdir6)
colnames(testdir6)[2] <- "InvExprVal"
colnames(testdir6)[3] <- "NatExprVal"
colnames(testdir6)[1] <- "Contig"
testdir6$InvUp <- testdir6$InvExprVal > testdir6$NatExprVal

sigOT0summary <- testdir6
head(sigOT0summary)
summary(sigOT0summary)
write.table(sigOT0summary, file="sigOrigin_popMeans_sumT0.txt", sep="\t")

####allowing for nearly equal things, sig origin####
#floating points in comparisons... don't do equals, just think of less than
close_enough <- function(x, y, tolerance=sqrt(.Machine$double.eps)) {
  abs(x - y) <= tolerance
}

sigOT0summary <- read.table(file="sigOrigin_popMeans_sumT0.txt", header=T)
sigOT0summaryNE <- sigOT0summary
sigOT0summaryNE$OriginsEq <- close_enough(sigOT0summaryNE$InvExprVal, sigOT0summaryNE$NatExprVal, tolerance=0.1)
sigOT0summaryNE$InvDefUp <- sigOT0summaryNE$InvUp 
sigOT0summaryNE[sigOT0summaryNE$OriginsEq==TRUE,]$InvDefUp <- FALSE
sigOT0summaryNE$InvDefDn <- sigOT0summaryNE$InvUp#FALSE means down
sigOT0summaryNE[sigOT0summaryNE$OriginsEq==TRUE,]$InvDefDn <- TRUE
write.table(sigOT0summaryNE, file="sigOrigin_popMeans_sumT0NE.txt", sep="\t")


####random effect pop####
PC1q_pop <- read.table("lme4_qval_PC1_pop.txt", header=T, sep="\t") 
head(PC1q_pop)
summary(PC1q_pop)
