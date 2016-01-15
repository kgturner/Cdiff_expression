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

##make tair list from my data
#split TAIRaccessions column
tair$TAIRaccessions <- as.character(tair$TAIRaccessions)
splitdat = do.call("rbind", strsplit(tair$TAIRaccessions, ","))
splitdat = data.frame(apply(splitdat, 2, as.character))
head(splitdat)
tair <- cbind(tair, splitdat)

#subset with significant range*treatment interaction or treatment
tair_dr <- subset(tair, intQ<=0.05 | trtQ<=0.05)
tair_dr <- droplevels(tair_dr)
tair_dr$X1 <- as.character(tair_dr$X1)
tair_dr$X2 <- as.character(tair_dr$X2)
tairList <- unique(c(tair_dr$X1, tair_dr$X2))

#check for duplicates
n_occur <- data.frame(table(tairList))
n_occur[n_occur$Freq > 1,]
tairList %in% n_occur$Var1[n_occur$Freq > 1]

##make Marchand list
Marchand$AGI_ID.from.HU <- as.character(Marchand$AGI_ID.from.HU)
Marchand$AGI_ID.from.HaT13l <- as.character(Marchand$AGI_ID.from.HaT13l)
MarchandList <- unique(c(Marchand$AGI_ID.from.HU, Marchand$AGI_ID.from.HaT13l))

#check for duplicates
n_occur <- data.frame(table(MarchandList))
n_occur[n_occur$Freq > 1,]
MarchandList %in% n_occur$Var1[n_occur$Freq > 1]

#compare lists
inGRNlist <- intersect(tairList, MarchandList) #52 unique TAIR accessions

##make table
#just accessions and descriptions that show up in both lists:
final <- subset(Marchand, AGI_ID.from.HU %in% inGRNlist |AGI_ID.from.HaT13l %in% inGRNlist)
final <- final[,4:9]
write.table(final, file="MarchandGRNoverlap.txt", sep="\t")

#checking duplication between columns - some easy editing in xcel
final.1 <- subset(final, AGI_ID.from.HU == AGI_ID.from.HaT13l) #some NAs in once col but not other
final.1 <- subset(final, name.from.Hu == name.from.HaT13l) #some NAs in once col but not other

