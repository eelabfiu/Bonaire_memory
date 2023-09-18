#Title: "Experimental set up, Environmental monitoring, Thermal stress assay"
#Author: "Serena Hackerott"
#Date: "08/29/2023"

#-------Set Up------------

####Load Packages####
library(viridis)
library(scales)
library(geodata)
library(ggplot2)

####Graphing Parameters####

####Colors by Groups

##Sites

##Genotypes

##Origins

##Colors by ordered factors
Site.colors.o<-c()

Gen.colors.o<-c()

Orig.colors.o<-c()


####Plot Feature Sizes
axis.title.sz=18
axis.txt.sz=14
leg.title.sz=15
leg.txt.sz=12
levels.sz=7
sig.sz=5
panel.lab.sz=25
point.sz=4
bar.sz=1
cap.sz=0.3


#-------Experimental Design------------

####Study Sites####


####Load Bonaire Map data
BES_remote <- "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_BES_1_sf.rds"

BES_rds <- file.path(tempdir(), "gadm36_BES_1_sf.rds")

if (toupper(Sys.info()["sysname"]) == "WINDOWS") {
  download.file(
    url = BES_remote,
    destfile = BES_rds,
    method = "wininet",
    mode = "wb"
  )
} else {
  download.file(
    url = BES_remote,
    destfile = BES_rds,
    method = "auto"
  )
}


##Read
BES.sf<-readRDS(BES_rds)

##Subset Bonaire
Bonaire.sf<-BES.sf[1,]


####Load Site Data
Sites<-read.csv("Data/BonaireSites.csv", header=TRUE)
str(Sites)

##Factors
Sites$Site<-factor(Sites$Site, levels=c("SS", "KL", "KR", "OL", "RP", "BN"), ordered=TRUE)
Sites$Type<-factor(Sites$Type, levels=c("Study", "Capital", "Harvest"), ordered=TRUE)


####Plot map of Bonaire with Study Sites
Study_Sites.plot<-ggplot(Bonaire.sf)+
  geom_sf(fill="white")+
  labs(x="Longitude", y="Latitude")+
  geom_point(data=Sites, 
             aes(x=Long, y=Lat, colour=Site, shape=Type), size=point.sz+3)+
  scale_colour_manual(values=c("#1F968BFF", "#32648EFF", "darkred", "grey25", "grey25", "grey25"))+
  scale_shape_manual(values=c(19, 17, 3))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz)); Study_Sites.plot

ggsave(filename="Figures/01_Design/Study_Sites.png", plot=Study_Sites.plot, dpi=300, width=8, height=12, units="in")




#-------Environmental Monitoring------------

####Load Data


#-------Acute Thermal Stress Assay------------

####Load Data
