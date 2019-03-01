library(popcycle)
library(dplyr)
library(lmodel2)

######################################
## DOWNLOAD THE RAW DATA using DAT ###
######################################
dat://fa2bb5981465ff583e091071f2f2c5e5f1b118c5219971a2dde107e4dec25e63

#Path to the raw data (DAT)
path.to.data <- "~/Documents/DATA/Codes/seaflow-virtualcore/seaflow-virtualcore-data/"



########################
### SEAFLOW ANALYSIS ###
########################

#Path to the Git repository
setwd("~/Documents/DATA/Codes/seaflow-virtualcore")



inst <- 740
print(inst)
list <- list.files(path=paste0(path.to.data,inst,"_caldata"),pattern=".evt$", recursive=T, full.names=T)
  if(inst == 740 | inst == 751) list <- list[-14]
  if(inst == 989)list <- list[-c(15:16)]




width <- 2500


ALL <- NULL

for (i in 1:length(list)){

#i <- 13
  print(list[i])

evt <- readSeaflow(list[i],transform=F)

  # Filtering particles detected by fsc_small
  evt. <- evt[evt$fsc_small > 1 | evt$D1 > 1 | evt$D2 > 1, ]
  evt. <- evt.[evt.$fsc_small > 1, ]

  aligned <- subset(evt., D2 < D1 + width & D1 < D2 + width )

  if(nrow(aligned)> 100000){aligned. <- sample_n(aligned, 100000)
  } else aligned. <- aligned


par(mfrow=c(1,1))
plot.cytogram(aligned., "fsc_small", "pe")
poly.beads <- getpoly(quiet=TRUE)
b <- subset(aligned,inout(aligned[,c("fsc_small", "pe")],poly=poly.beads, bound=TRUE, quiet=TRUE))
b <- subset(b, D1 < 60000 | D2 < 60000)

plot.cytogram(b, "fsc_small", "D1")
poly.opp <- getpoly(quiet=TRUE)
opp <- subset(b,inout(b[,c("fsc_small", "D1")],poly=poly.opp, bound=TRUE, quiet=TRUE))


png(paste0("1.bead_calibration/",inst,"_caldata/",basename(list[i]),".png"),width=9, height=12, unit='in', res=100)
  par(mfrow=c(3,2))
  plot.cytogram(aligned., "fsc_small", "pe"); polygon(poly.beads, border=2, lwd=2)
  plot.cytogram(aligned., "fsc_small", "chl_small")
  plot.cytogram(b, "fsc_small", "D1"); polygon(poly.opp, border=2, lwd=2)
  plot.cytogram(b, "fsc_small", "D2")
  plot.cytogram(opp, "fsc_small", "pe")
  plot.cytogram(opp, "fsc_small", "chl_small")
dev.off()

beads <- nrow(opp)
fsc.med <- quantile(opp$fsc_small, 0.5)
fsc.low <- quantile(opp$fsc_small, 0.025)
fsc.high <- quantile(opp$fsc_small, 0.975)
d1.med <- quantile(opp$D1, 0.5)
d1.low <- quantile(opp$D1, 0.025)
d1.high <- quantile(opp$D1, 0.975)
d2.med <- quantile(opp$D2, 0.5)
d2.low <- quantile(opp$D2, 0.025)
d2.high <- quantile(opp$D2, 0.975)

all <- data.frame(cbind(file=as.character(list[i]), beads, fsc.med, fsc.low, fsc.high, d1.med, d1.low, d1.high, d2.med, d2.high, d2.low))

ALL <- rbind(ALL, all)

}

write.csv(ALL, paste0("1.bead_calibration/",inst,"-summary.csv"), quote=F, row.names=F)




############
## FIGURE ##
############
library(popcycle)
library(tidyverse)
path.to.data <- "~/Documents/DATA/Codes/seaflow-virtualcore/seaflow-virtualcore-data/"
inst <- 751
list <- list.files(path=paste0(path.to.data,inst,"_caldata"),pattern=".evt$", recursive=T, full.names=T)

size <- c(0.3,0.5,0.7,1,1.8,3.1,5.7)
id <- c(7,8,9,11,10,12,13) #oligo
#id <- c(15,1,2,3,16,4,5) #coastal

BEADS <- NULL
for(i in 1:length(id)){
  evt <- readSeaflow(list[id[i]],transform=T)
  evt. <- subset(evt, fsc_small > 1)
  evt. <- subset(evt., D1 < max(evt.$D1) & D2 < max(evt.$D2))
  plot.cytogram(evt., para.x='fsc_small', para.y='pe')
  poly.beads <- getpoly(quiet=TRUE)
  beads <- subset(evt.,inout(evt.[,c("fsc_small", "pe")],poly=poly.beads, bound=TRUE, quiet=TRUE))
  beads$size <- size[i]
  BEADS <- rbind(BEADS, beads)
}

write.csv(BEADS, paste0("1.bead_calibration/",inst,"-EVT.csv"), quote=F, row.names=F)

p <- plot_cytogram(BEADS, para.x='fsc_small', para.y='D1', bins=200)
p + geom_point(aes(x=c(10,100,500,750,1000,1200), y=c(3,9,50,70,100,300)))
