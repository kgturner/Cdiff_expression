#C.diffusa expression - TAIR tables of sig genes
#3/24/15

TR001 <- read.table(file = "~/GitHub/Cdiff_expression/GOanalysis/unique.sorted.out.ath_Cendif1.unigenes_GO", sep=" ") #go annots from BLAST
US022 <- read.table(file = "~/GitHub/Cdiff_expression/GOanalysis/unique.sorted.out.ath_Cendif2.unigenes_GO", sep=" ")
GOboth <- read.table("~/GitHub/Cdiff_expression/GOanalysis/GOmap_Cendif1_Cendif2.txt_flip", sep="\t")

PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") #contigs and q values


table <- merge(PC1q, GOboth, by.x="Contig", by.y="V1", all.x=TRUE)

TR001a <- TR001[,1:2]
TR001a <- unique(TR001a)
summary(TR001a)

US022a <- US022[,1:2]
US022a <- unique(US022a)

bothA <- merge(TR001a, US022a, all=TRUE)
# bothA <- unique(bothA)

table2 <- merge(table, bothA, by.x="Contig", by.y="V2", all.x=TRUE )
head(table2)

table3 <- subset(table2, intQsig==TRUE|covQsig==TRUE|originQsig==TRUE|trtQsig==TRUE, select=c(Contig, V1, intQ, originQ, trtQ, covQ, V2))
#should be fewer rows than 9578+1111+585+227???
summary(table3)
subset(table3, Contig=="Contig100")
#ah, multiple TAIR accessions for each Contig...

# library(reshape2)
# table4 <- reshape(table3, direction="wide", idvar="Contig", v.names="V1", timevar=rownames(table3))

# library(plyr)
# ## Add a medication index
# table3_index <- ddply(table3, .(Contig), mutate, 
#                          index = paste0('AT', 1:length(Contig)))    
# table4 <- dcast(table3_index, Contig~ index, value.var = 'TAIRacc')

library(splitstackshape)
# test <- head(table3)
# table3_id <- getanID(table3, "Contig")
table3_id <- dcast.data.table(getanID(table3, "Contig"), Contig ~ .id, value.var = "V1")
colnames(table3_id) <- c("Contig", "AT1", "AT2")
table3_id$TAIRaccessions <- paste0(table3_id$AT1,",",table3_id$AT2)

table4 <- merge(table3, table3_id, by="Contig", all.y=TRUE)
table4 <- subset(table4, select=c(Contig, TAIRaccessions, intQ, originQ, trtQ, covQ, V2))
table5 <- unique(table4)
table5$TAIRaccessions <- as.factor(table5$TAIRaccessions)
colnames(table5)[7] <- "GOterms"
rownames(table5) <- table5$Contig

write.table(table5, file="LRTsandTAIRaccessions.txt", sep="\t")
