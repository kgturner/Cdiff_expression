#C.diffusa expression - pick out differentially expressed genes
#7/20/15

#load data
lrtair <- read.table("LRTsandTAIRaccessions.txt", header=T, sep="\t") #lrt pvals and tair accessions

#up and down data

####constitutive differences, timepoint 0 control####
sigOT0summaryNE <- read.table(file="sigOrigin_popMeans_sumT0NE.txt", header=T)

#merge, add fold differences
sigOT0tair <- merge(sigOT0summaryNE, lrtair, all.x=TRUE)
sigOT0tair$folddiff <- sigOT0tair$InvExprVal/sigOT0tair$NatExprVal #fold differences

#log-fold difference
extremefold <- subset(sigOT0tair, (folddiff>=1.5|folddiff<=0.67)&TAIRaccessions!="NA,NA") 
#extreme fold difference >1.5 (1.5 fold increase in inv rel to nat) or <.67 (1.5 fold decrease in inv rel to nat);  
#only ones with TAIR accessions:
# AT3G16630 >1.5 fold higher in inv
# AT5G67280 >1.5 fold higher in inv
# AT2G46100 >1.5 fold lower in inv

#sig genes in particular pathways
#etheylene, down in inv
# GO:0007165  signal transduction
GO:0009873  ethylene-activated signaling pathway
GO:0071369	cellular response to ethylene stimulus

eth <- subset(sigOT0tair, InvDefDn=T&TAIRaccessions!="NA,NA") 
eth2 <- grep("0009873|0071369", eth$GOterms, perl=TRUE, value=FALSE)
ethsub <- subset(eth, row.names(eth)%in%eth2)


str(eth$GOterms)

# ethsub <- for (i in 1:nrow(eth)) {
#   if (grepl("0007165|0009873|0071369",eth[i,"GOterms"])){
#   return(eth[i,])}
#   
# }







grep("a+", c("abc", "def", "cba a", "aa"), perl=TRUE, value=TRUE)





####induced difference, timepoint 2, drought####
int2dsummaryNE <- read.table(file="sigint_popMeans_sumT2droughtNE.txt", header=T)

#merge, add fold differences
sigintT2tair <- merge(int2dsummaryNE, lrtair, all.x=TRUE)
sigintT2tair$folddiff <- sigintT2tair$InvExprVal/sigintT2tair$NatExprVal #fold differences

extremefoldT2 <- subset(sigintT2tair, (folddiff>=1.2|folddiff<=0.84)&TAIRaccessions!="NA,NA") 
#extreme fold difference >1.5 (1.5 fold increase in inv rel to nat) or <.67 (1.5 fold decrease in inv rel to nat); 
#only ones with TAIR accessions:
#none at 1.5, change to 1.2, only one contig
# AT4G30410,AT4G23420 >1.2 fold increase in inv rel to nat

#
