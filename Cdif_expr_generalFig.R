#generalized theorethical figure for Cdif gene expression paper
#2/9/2016

library(plyr)
library(ggplot2)

dat<- read.table("BSfig_data.txt", header=T, sep="\t",quote='"', row.names=1) #"BSfig_data.txt"

####generate data####
# n <- 20 #how many numbers you want
# 
# m <- 1 #mean
# s <- 1 #std dev
# 
# L <- .2 #expected min
# U <- .8 #expected max
# 
# p_L <- pnorm(L, mean=m, sd=s)
# p_U <- pnorm(U, mean=m, sd=s) 
# x <- qnorm(runif(n, p_L, p_U), mean=m, sd=s) 


#or
rnorm(n=20, m=1, sd=0.25) 

####data vectors (for 40 pops)####
# tol_nochange <- c(rnorm(n=20, m=2, sd=0.25),rnorm(n=20, m=2, sd=0.25))
# size_nochange <- c(rnorm(n=20, m=2, sd=0.25),rnorm(n=20, m=2, sd=0.25))

size_decr <- c(rnorm(n=20, m=2, sd=0.25),rnorm(n=20, m=1, sd=0.25))
# size_incr <- c(rnorm(n=20, m=2, sd=0.25),rnorm(n=20, m=3, sd=0.25))

tol_incr <- c(rnorm(n=20, m=2, sd=0.25),rnorm(n=20, m=3, sd=0.25))
# tol_decr <- c(rnorm(n=20, m=2, sd=0.25),rnorm(n=20, m=1, sd=0.25))

####unidirectional trade off for incr tolerance in inv####
dat_incrTol <- subset(dat, Group%in%c("Native", "Invasive A"))

dat_incrTol$ControlSize <- size_decr
dat_incrTol$StressTolerance <- tol_incr

colnames(dat_incrTol)[1] <- "Range"

pTradeTol <- ggplot(dat_incrTol,aes(x=StressTolerance, y=ControlSize)) +  
  geom_point(aes(shape=Range, color=Range), size=4)+
  geom_smooth(method=glm, se=TRUE,color="black", size=2 )+
  #   coord_cartesian(ylim = c(0, 4), xlim= c(0, 4)) +
  theme_bw()+
  theme(axis.ticks = element_blank(), axis.text = element_blank())+
  xlab("Breadth of environmental tolerance") + 
  ylab("Magnitude of expression response to drought") +
#   annotate(geom="text", x=1.5, y=2.4, label="(b)",fontface="bold", size=3)+
  theme(legend.justification=c(0.95,0.95), legend.position=c(1,1),
                legend.title = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12))

pTradeTol

##curves?
pCurve <- ggplot(data = dat_incrTol, mapping = aes(x=StressTolerance, y=ControlSize)) +
  layer(geom = "line") +
  layer(geom = "area", #mapping = aes(x = ifelse(StressTolerance>-1.2 & x<1.1 , x, 0)
        geom_params = list(fill = "red", alpha = 0.5)) +
  scale_y_continuous(limits = c(0, max(dat_incrTol$ControlSize)))
pCurve

####figure####
png("Cdifexpr_Fig4.png",width=400, height = 400, pointsize = 16)
# pdf("Cdifexpr_Fig4.pdf", useDingbats=FALSE, width=4.4, height=4.4, pointsize = 16) #MolEcol sizes 3.149, 4.4 or 6.65

pTradeTol
dev.off()



# ####multipanel plot####
# # pdf("InvaGenTradeoffFig.pdf", useDingbats=FALSE, width=6.65, height=9, pointsize = 12) #MolEcol sizes 3.149, 4.4 or 6.65
# # multiplot(pSizeB, pTolB, pTrade, cols=2)
# # # legend("top", c("Invasive C. diffusa","Native C. diffusa", "Native C. stoebe"), 
# # #        pch=c(16,17,15), fill=origincol,  bg="white", title = "Sampled populations", cex=1)
# # dev.off()

# ####multiplot####
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   require(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#     
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }