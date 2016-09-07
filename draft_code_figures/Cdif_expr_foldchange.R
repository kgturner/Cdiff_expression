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
#for InvDefDn, FALSE means down
# GO:0009873  ethylene-activated signaling pathway
# GO:0071369	cellular response to ethylene stimulus

eth <- subset(sigOT0tair, InvDefDn==FALSE&TAIRaccessions!="NA,NA") 
eth2 <- eth[grep("0009873|0071369", eth$GOterms, perl=TRUE, value=FALSE),]

# AT3G46060,AT5G19440
# AT1G64060

# str(eth$GOterms)
# 
# # ethsub <- for (i in 1:nrow(eth)) {
# #   if (grepl("0007165|0009873|0071369",eth[i,"GOterms"])){
# #   return(eth[i,])}
# #   
# # }
# 
# grep("a+", c("abc", "def", "cba a", "aa"), perl=TRUE, value=TRUE)

#cp, up in inv
# GO:0030095  chloroplast photosystem II
# GO:0009654	photosystem II oxygen evolving complex
# GO:0009521	photosystem
cp <- subset(sigOT0tair, InvDefUp==TRUE&TAIRaccessions!="NA,NA") 
cp2 <- cp[grep("0030095|0009654|0009521", cp$GOterms, perl=TRUE, value=FALSE),]

# AT1G10360
# AT3G07890
# AT4G15510

#mt, up in inv
# GO:0007005  mitochondrion organization
mt <- subset(sigOT0tair, InvDefUp==TRUE&TAIRaccessions!="NA,NA") 
mt2 <- mt[grep("0007005", mt$GOterms, perl=TRUE, value=FALSE),]

# AT1G17020,AT1G14830

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
 
#sig genes in particular pathways
#mt, up in inv
# GO:0045039  protein import into mitochondrial inner membrane
# GO:0090151  establishment of protein localization to mitochondrial membrane
# GO:0007006  mitochondrial membrane organization
# GO:0006626  protein targeting to mitochondrion
# GO:0070585  protein localization to mitochondrion
# GO:0072655  establishment of protein localization to mitochondrion
# GO:0005740  mitochondrial envelope
# GO:0005758  mitochondrial intermembrane space

mtstar <- subset(sigintT2tair, InvDefUp==TRUE&TAIRaccessions!="NA,NA") 
mtstar2 <- mtstar[grep("0045039|0090151|0007006|0006626|0070585|0072655|0005740|0005758", mtstar$GOterms, perl=TRUE, value=FALSE),]

# AT5G42300,AT3G15640
# AT5G50810
# AT2G29530
# AT1G80230
# AT1G52710

#glycolipids, down in inv
#for InvDefDn, FALSE means down
# GO:0006664  glycolipid metabolic process
# GO:0009247  glycolipid biosynthetic process


gly <- subset(sigintT2tair, InvDefDn==FALSE&TAIRaccessions!="NA,NA") 
gly2 <- gly[grep("0006664|0009247", gly$GOterms, perl=TRUE, value=FALSE),]

# AT4G33030

#vacuole, down in inv
#for InvDefDn, FALSE means down
# GO:0005774  vacuolar membrane
# GO:0005773  vacuole
# GO:0044437  vacuolar part

vac <- subset(sigintT2tair, InvDefDn==FALSE&TAIRaccessions!="NA,NA") 
vac2 <- vac[grep("0005774|0005773|0044437", vac$GOterms, perl=TRUE, value=FALSE),]

# AT1G04640,AT2G31610
# AT2G07050
# AT4G39660
# AT3G28480

