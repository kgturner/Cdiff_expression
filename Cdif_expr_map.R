#Cdif gene expression - population map
#2/3/2016


library(rgdal) # Commands for reprojecting the vector data.
library(rworldmap) # Recently updated mapping program.
library(rworldxtra) # Add-ons for rworldmap.

####set up map####
projectionCRS <- CRS("+proj=laea +lon_0=0.001 +lat_0=89.999 +ellps=sphere") #the ellps 'sphere' has a radius of 6370997.0m
par(mai=c(0,0,0.2,0)) #,xaxs="i",yaxs="i"
sPDF <- getMap()[-which(getMap()$ADMIN=='Antarctica')] 
sPDF <- spTransform(sPDF, CRS=projectionCRS)
setLims <- TRUE #FALSE back to whole world
#setLims <- FALSE
if ( !setLims )
{
  xlim <- ylim <- NA
} else
{
  ### TRY FIDDLING WITH THESE LIMITS ###
  xlimUnproj <- c(-52,120)
  ylimUnproj <- c(10,30)
  sPointsLims <- data.frame(x=xlimUnproj, y=ylimUnproj)
  coordinates(sPointsLims) = c("x", "y")
  proj4string(sPointsLims) <- CRS("+proj=longlat +ellps=WGS84")
  sPointsLims <- spTransform(sPointsLims, CRS=projectionCRS)
  xlim <- coordinates(sPointsLims)[,"x"]
  ylim <- coordinates(sPointsLims)[,"y"]  
}

# sPDF <- getMap()
# #list of country names
# sPDF$ADMIN
#setup a color code column filled with numbers
sPDF$colCode <- 4

#set codes for specified countries
sPDF$colCode[ which(sPDF$ADMIN %in% c("Canada","United States of America"))] <- 1
sPDF$colCode[ which(sPDF$ADMIN %in% c("Armenia","Azerbaijan", "Bulgaria", "Georgia", 
                                      "Greece", "Moldova", "Romania","Russia", "Turkey",
                                      "Ukraine", "Serbia"))] <- 2
sPDF$colCode[ which(sPDF$ADMIN %in% c("Poland", "Belarus", "Italy", "Syria", "Czech Republic",
                                      "Estonia", "Switzerland","Latvia","Lithuania", 
                                      "Slovenia", "Serbia","Austria","Belgium", "France",
                                      "Germany","Hungary","Luxembourg","Norway","Slovakia",
                                      "Spain", "United Kingdom", "Kazakhstan", "Turkmenistan", "China"))] <- 3

#create a colour palette - note for each value not for each country
colourPalette <- c("#F8766D","#00BFC4", "cadetblue","lightgray") #inv, nat, present/naturalized, extra countries

# spName <- plotmath(italic("Centaurea diffusa"))

#points
pop <- read.table("Popcoord.txt", header=TRUE, stringsAsFactor=FALSE)
pop <- subset(pop, Pop%in%c("BG001", "RU008", "TR001", "CA001", "US001", "US003" ))

pop$Origin <- "Native"
pop[pop$Pop%in%c("CA001", "US001", "US003"),]$Origin <- "Invasive"

pop$Pop <- as.factor(pop$Pop)
pop$Origin <- as.factor(pop$Origin)
pop$Latitude <- as.numeric(pop$Latitude)
pop$Longitude <- as.numeric(pop$Longitude)

pop$pch <- 1 #for invasives
pop[pop$Origin %in% "Native",]$pch <- 17

coordinates(pop) = c("Longitude", "Latitude")
proj4string(pop) <- CRS("+proj=longlat +ellps=WGS84")
sPointsDF <- spTransform(pop, CRS=projectionCRS)

#lat markings...
markings <- data.frame(Latitude=as.numeric(c(75,60,45,30,15,85,85)), Longitude=as.numeric(c(-45,-45,-45,-45,-45,0,180)),name=c("75", "60","45","30","15","0","180"))
coordinates(markings) = c("Longitude", "Latitude")
proj4string(markings) <- CRS("+proj=longlat +ellps=WGS84")
sPointsDFmark <- spTransform(markings, CRS=projectionCRS)

####plot map####
# pdf("KTurnerFig1.pdf", useDingbats=FALSE, width=6.65, height = 5, pointsize = 12) #4.4 or 6.65
png("Cdif_exprMap.png", width=665, height = 500, pointsize = 12)
# svg("KTurnerFig1.svg", width=6.65, height = 5, pointsize = 12)

par(mar=c(0,0,0,0))
mapCountryData(sPDF, nameColumnToPlot="colCode", mapTitle=NA,
               colourPalette=colourPalette, borderCol ='gray24', addLegend = FALSE,
               xlim=xlim, ylim=ylim, catMethod=c(0,1,2,3,4))
#note that catMethod defines the breaks and values go in a category if they are <= upper end
#mapTitle=bquote(Global~range~of~italic(Centaurea)~italic(diffusa)) 
points(sPointsDF, pch=pop$pch, cex=2, lwd=2)

llgridlines(sPDF, easts=c(-90,-180,0,90,180), norths=seq(0,90,by=15), 
            plotLabels=FALSE, ndiscr=1000) #ndiscr=num points in lines
text(sPointsDFmark, labels = sPointsDFmark$name, cex=1) #pch2 for triangles

legend("topright", c("Invasive C. diffusa","Native C. diffusa"), 
       pch=c(1,17),  bg="white", title = "Sampled populations", cex=1)
legend("bottomleft", c("Invasive", "Native","Naturalized"), fill=colourPalette,
       title="Ranges of C. diffusa", bg="white", cex=1)
box(lty="solid", col = "black")

dev.off()