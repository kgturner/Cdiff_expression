#C.diffusa expression - time course
#2/16/15
#based on Gillespie et al., 2010. Github repo at https://github.com/csgillespie/bmc-microarray

#can't handle complex model structure, with Origin and Trt????

source("http://bioconductor.org/biocLite.R")
biocLite(c("timecourse"))

library(timecourse)

#load data
PC1q_intsigdf <- read.table("PC1_sigint_df.txt", header=T)
PC1q_Osigdf <- read.table("PC1_sigOrigin_df.txt", header=T)
PC1q_sigdf <- read.table("PC1_sigPC1_df.txt", header=T)
PC1q_trtsigdf <- read.table("PC1_sigtrt_df.txt", header=T)

#order by condition rep, and time
tc_intdf <- PC1q_intsigdf
tc_intdf$OTrt <- as.factor(paste0(tc_intdf$Origin, "_", tc_intdf$Trt))

temp <- strsplit(as.character(tc_intdf$TrtPool), "")
temp2 <- matrix(unlist(temp), ncol=2, byrow=TRUE)
tc_intdf$Pool <- as.factor(temp2[,2])

tc_intdf <- tc_intdf[order(tc_intdf$OTrt, tc_intdf$TrtPool, tc_intdf$Tmpt), ]
tc_intdf[,c(1:14, 242:243)]

#calculate the Hotelling statistic T˜ 2 via
c.grp = as.character(tc_intdf$OTrt) #condition
t.grp = as.numeric(tc_intdf$Tmpt) #time
r.grp = as.character(tc_intdf$TrtPool) #replicate

tc_intM <- as.matrix(t(tc_intdf[,c(15:241)]))
size = matrix(6, nrow = 227, ncol = 4) #36 is number of levels from PopTrtPool?, nrow is number of genes in subsetted dataset, 2 is number of treatments
# size <- rep(3, 2000)

MB.multi = mb.MANOVA(tc_intM, times = 3, D=4,
                size = size, condition.grp = c.grp)
summary(MB.multi)

# The top (say) one hundred genes can be extracted via
gene_positions = MB.multi$pos.HotellingT2
gnames = rownames(tc_intM)
gene_probes = gnames[gene_positions]
# The expression profiles can also be easily obtained. The profile for the top ranked expression is
# found using
plotProfile(MB.multi)
