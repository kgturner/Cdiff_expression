#C. diffusa expression - comparing to results from other papers
#1/14/2016

####compare GO terms to InvSyn results####
#Hodgins et al., Mol Ecol, 2015
#to see how lists produced, see Cdif_expr_topgo.R
#GO terms enriched for rapidly evolving genes in invasive diffusa
invsyn <- read.table("GOterm_InvSyn.txt", header=T, sep="\t")

intpcRes$GO.ID %in% invsyn$GO.ID
opcRes$GO.ID %in% invsyn$GO.ID
invUp2d_pcRes$GO.ID %in% invsyn$GO.ID
invDn2d_pcRes$GO.ID %in% invsyn$GO.ID
invUpoT0_pcRes$GO.ID %in% invsyn$GO.ID
invDnoT0_pcRes$GO.ID %in% invsyn$GO.ID

PC1pcRes <- read.table(file="GOresults_sigPC1_pc.txt", header=T)
PC1pcRes$GO.ID %in% invsyn$GO.ID

MFpcRes$GO.ID %in% invsyn$GO.ID
CCintpcRes$GO.ID %in% invsyn$GO.ID
MFopcRes$GO.ID %in% invsyn$GO.ID
CCopcRes$GO.ID %in% invsyn$GO.ID

####Compare genes to sunflower gene regulatory network for drought####
#Marchand et al., New Phyt, 2014

#load data
tair <- read.table("LRTsandTAIRaccessions.txt", header=T, sep="\t")
Marchand <- read.table ("nph12818TableS3.txt", header=T, sep="\t")

#split TAIRaccessions column
tair$TAIRaccessions <- as.character(tair$TAIRaccessions)
# unlist(strsplit("a.b.c", ".", fixed = TRUE))

splitdat = do.call("rbind", strsplit(tair$TAIRaccessions, ","))
splitdat = data.frame(apply(splitdat, 2, as.character))
# names(splitdat) = paste("trial", 1:4, sep = "")
head(splitdat)

tair <- cbind(tair, splitdat)
#subset with significant range*treatment interaction or treatment
tair_dr <- subset(tair, intQ<=0.05 | trtQ<=0.05)

#compare lists
#AT codes from Marchand in $AGI_ID.from.HU and $AGI_ID.from.HaT13l
#and in tair from $X1 and $X2

inGRN <- subset(tair_dr, X1 %in% Marchand$AGI_ID.from.HU)
inGRN <- rbind(inGRN, subset(tair, X2 %in% Marchand$AGI_ID.from.HU))
inGRN <- rbind(inGRN, subset(tair, X1 %in% Marchand$AGI_ID.from.HaT13l))
inGRN <- rbind(inGRN, subset(tair, X2 %in% Marchand$AGI_ID.from.HaT13l))
inGRN <- droplevels(inGRN)

inGRNlist <- unique(c(levels(inGRN$X1),levels(inGRN$X2)))
# inGRNlist <- unique(inGRNlist)
inGRNlist <- inGRNlist[1:89] #remove NA

MarchandOverlap <- subset(Marchand, AGI_ID.from.HU %in% inGRNlist | AGI_ID.from.HaT13l %in% inGRNlist)
MarchandOverlap <- MarchandOverlap[,1:9]
tairOverlap <- subset(tair_dr, X1 %in% inGRNlist | X2 %in% inGRNlist)
tairOverlap <- droplevels(tairOverlap)

#just accessions and descriptions that show up in both lists:
#not working yet...
MarchandComp <- as.data.frame(inGRNlist)
final <- merge(MarchandComp, Marchand, by.x="inGRNlist", by.y="AGI_ID.from.HU", all.x=T)
final <- merge(final, Marchand, by.x="inGRNlist", by.y="AGI_ID.from.HaT13l", all.x=T)
# final$name.from.Hu %in% final$name.from.HaT13l


#attepmts to make fancier table
colnames(tairOverlap)[8] <- "AT1"
colnames(tairOverlap)[9] <- "AT2"
tairMarch.1 <- tairOverlap[,c(1,8)]
tairMarch.2<- tairOverlap[,c(1,9)]
tairMarch.2 <- subset(tairMarch.2, AT2!="NA")
colnames(tairMarch.1)[2] <- "TAIRaccession"
colnames(tairMarch.2)[2] <- "TAIRaccession"
tairMarch <- unique(rbind(tairMarch.1, tairMarch.2))
rownames(tairMarch) <- tairMarch$TAIRaccession


March <- MarchandOverlap[,c(4:8)]


tairMarch <- merge(tairMarch, March, by.x="AT1", by.y="AGI_ID.from.HU", all.x=T)
tairMarch2 <- merge(tairMarch, March, by.x="AT2", by.y="AGI_ID.from.HU", all.x=T)
# tairMarch3<- merge(tairMarch, March, by.x="AT1", by.y="AGI_ID.from.HaT13l", all.x=T)
# tairMarch4<- merge(tairMarch, March, by.x="AT2", by.y="AGI_ID.from.HaT13l", all.x=T)

tairMarchAll <- merge(tairMarch, tairMarch2,by="Contig")# 
# # tairMarchall <- merge(tairMarch, tairMarch2, by="Contig", all=T)
# tairMarchall <- 

squish <- function(dat=dat, cols=cols, newcol="newcol"){
  
  dat$temp <- NA
  
  for (i in dat[,cols]){
    ## where is temp NA?
    ss <- which(is.na(dat$temp))
    ## what are the values of i there? put them in temp
    dat$temp[ss] <- i[ss] 
  }
  
  names(dat)[names(dat)=="temp"] <- newcol
  
  return(subset(dat,select=ncol(dat)))
}
test <- squish(dat=tairMarch2, cols=c(4,8), newcol="alt")
