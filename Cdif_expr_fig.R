#C.diffusa expression - figures
#2/26/15



#load data
PC1q <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") #contigs and q values
PC1p <- read.table("lme4_PC1_pvalues.txt", header=T, sep="\t") #contigs and p values
PC1q_intsigdf <- read.table("PC1_sigint_df.txt", header=T)
PC1q_Osigdf <- read.table("PC1_sigOrigin_df.txt", header=T)
PC1q_sigdf <- read.table("PC1_sigPC1_df.txt", header=T)
PC1q_trtsigdf <- read.table("PC1_sigtrt_df.txt", header=T)

####venn diagram####
#http://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html
install.packages('VennDiagram')
library(VennDiagram)

