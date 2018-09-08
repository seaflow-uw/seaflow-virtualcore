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


png(paste0("/1.bead_calibration/",inst,"_caldata/",basename(list[i]),".png"),width=9, height=12, unit='in', res=100)
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

write.csv(ALL, paste0("/1.bead_calibration/",inst,"-summary.csv"), quote=F, row.names=F)








##############
### SLOPES ###
##############
SLOPES <- NULL

png(paste0("/1.bead_calibration/ALL_MERGED_Slopes.png"),width=12, height=15, unit='in', res=100)
par(mfrow=c(3,2))

  for(ins in c(989,751,740)){
  print(ins)

  ALL <- read.csv(paste0("/1.bead_calibration/",ins,"-summary.csv"))
    ALL$pmt <- 'coastal'
    id <- which((grepl("oligo", ALL$file))==T)
    ALL[id, "pmt"] <- 'oligo'
    ALL$size <- as.numeric(sub("um.evt","",t(data.frame(list(strsplit(as.character(ALL$file),"_"))))[,5], ))
    ALL <- ALL[order(ALL$size),]


     df.o <- subset(ALL, pmt == 'oligo'  & size > 0.18)
     ref.o <- subset(df.o, size ==1.00)
     df.c <- subset(ALL, pmt == 'coastal'  & size > 0.18)
     ref.c <- subset(df.c, size ==1.00)

     if(ins == 740){
          comp.fsc <- mean(df.o$fsc.med[-1] - df.c$fsc.med)
          comp.d1 <- mean(df.o$d1.med[-1] - df.c$d1.med)
          comp.d2 <- mean(df.o$d2.med[-1] - df.c$d2.med)
          }else{
         comp.fsc <- mean(df.o$fsc.med - df.c$fsc.med)
         comp.d1 <- mean(df.o$d1.med - df.c$d1.med)
         comp.d2 <- mean(df.o$d2.med - df.c$d2.med)
          }

      DF <- data.frame(cbind(fsc=c(df.o$fsc.med, df.c$fsc.med+comp.fsc),#df.o$fsc.low,df.o$fsc.high,df.c$fsc.low+comp.fsc,df.c$fsc.high+comp.fsc),
                              d1=c(df.o$d1.med, df.c$d1.med+comp.d1),#df.o$d1.low,df.o$d1.high,df.c$d1.low+comp.d1,df.c$d1.high+comp.d1),
                              d2=c(df.o$d2.med, df.c$d2.med+comp.d2),#df.o$d2.low,df.o$d2.high,df.c$d2.low+comp.d2,df.c$d2.high+comp.d2),
                              size=c(df.o$size,df.c$size)))

      # fit <- lm(d1 ~ stats:::poly(fsc,2, raw=T),data=DF)
      #   summary(fit)
      # pred <- predict(fit, data.frame(fsc=seq(0,65000, by=1000)))
      # plot(DF$fsc, DF$d1, xlim=c(0,65000), ylim=c(0,65000),main=paste(ins))
      # lines(seq(0,65000, by=1000),pred, lwd=2, col=2)

  #SMALL
  sm <- subset(DF, size <= 1)
  sm.lm.d1 <- lmodel2(d1 ~ fsc, data=sm, range.x="interval", range.y="interval", nperm=99)
  sm.lm.d2 <- lmodel2(d2 ~ fsc, data=sm, range.x="interval", range.y="interval", nperm=99)

  # #1 MICRON BEADS
  # ref <- subset(DF, size == 1)
  # ref.fsc <- mean(ref$fsc)
  # #BREAK
  # br <- subset(DF, size == 1.83 | size ==3.10)
  # break.fsc <- mean(br$fsc)

  #LARGE
  lg <- subset(DF, size >= 1)
  lg.lm.d1 <- lmodel2(d1 ~ fsc, data=lg, range.x="relative", range.y="relative", nperm=99)
  lg.lm.d2 <- lmodel2(d2 ~ fsc, data=lg, range.x="relative", range.y="relative", nperm=99)

  notch.small.D1 <-  round(sm.lm.d1$regression.results$Slope[1],3)
  notch.small.D2 <-  round(sm.lm.d2$regression.results$Slope[1],3)
  notch.large.D1 <-  round(lg.lm.d1$regression.results$Slope[1],3)
  notch.large.D2 <-  round(lg.lm.d2$regression.results$Slope[1],3)
  intersect.small.D1 <- round(sm.lm.d1$regression.results$Intercept[1])
  intersect.small.D2 <- round(sm.lm.d2$regression.results$Intercept[1])
  intersect.large.D1 <- round(lg.lm.d1$regression.results$Intercept[1])
  intersect.large.D2 <- round(lg.lm.d2$regression.results$Intercept[1])

  notch.small.D1_97.5 <-  round(sm.lm.d1$confidence.intervals[1,5],3)
  notch.small.D2_97.5 <-  round(sm.lm.d2$confidence.intervals[1,5],3)
  notch.large.D1_97.5 <-  round(lg.lm.d1$confidence.intervals[1,5],3)
  notch.large.D2_97.5 <-  round(lg.lm.d2$confidence.intervals[1,5],3)
  intersect.small.D1_97.5 <- round(sm.lm.d1$confidence.intervals[1,2])
  intersect.small.D2_97.5 <- round(sm.lm.d2$confidence.intervals[1,2])
  intersect.large.D1_97.5 <- round(lg.lm.d1$confidence.intervals[1,2])
  intersect.large.D2_97.5 <- round(lg.lm.d2$confidence.intervals[1,2])

  notch.small.D1_2.5 <-  round(sm.lm.d1$confidence.intervals[1,4],3)
  notch.small.D2_2.5 <-  round(sm.lm.d2$confidence.intervals[1,4],3)
  notch.large.D1_2.5 <-  round(lg.lm.d1$confidence.intervals[1,4],3)
  notch.large.D2_2.5 <-  round(lg.lm.d2$confidence.intervals[1,4],3)
  intersect.small.D1_2.5 <- round(sm.lm.d1$confidence.intervals[1,3])
  intersect.small.D2_2.5 <- round(sm.lm.d2$confidence.intervals[1,3])
  intersect.large.D1_2.5 <- round(lg.lm.d1$confidence.intervals[1,3])
  intersect.large.D2_2.5 <- round(lg.lm.d2$confidence.intervals[1,3])

  slope <- data.frame(cbind(ins, notch.small.D1,notch.small.D2,notch.large.D1,notch.large.D2,
                                intersect.small.D1,intersect.small.D2,intersect.large.D1,intersect.large.D2,
                            notch.small.D1_97.5,notch.small.D2_97.5,notch.large.D1_97.5,notch.large.D2_97.5,
                                intersect.small.D1_97.5,intersect.small.D2_97.5,intersect.large.D1_97.5,intersect.large.D2_97.5,
                            notch.small.D1_2.5,notch.small.D2_2.5,notch.large.D1_2.5,notch.large.D2_2.5,
                              intersect.small.D1_2.5,intersect.small.D2_2.5,intersect.large.D1_2.5,intersect.large.D2_2.5))
  SLOPES <- rbind(SLOPES, slope)

    plot(DF$fsc, DF$d1, xlim=c(0,65000), ylim=c(0,65000),main=paste(ins))
      points(df.c$fsc.med+comp.fsc,df.c$d1.med+comp.d1, col=2)
      points(ref.o$fsc.med, ref.o$d1.med, pch=16)
      points(ref.c$fsc.med+comp.fsc, ref.c$d1.med+comp.d1, col=2,pch=16)
      abline(b=notch.small.D1, a=intersect.small.D1,col=4)
      abline(b=notch.large.D1, a=intersect.large.D1,col=3)
        abline(b=notch.small.D1_97.5, a=intersect.small.D1_97.5,col='grey')
        abline(b=notch.small.D1_2.5, a=intersect.small.D1_2.5,col='grey')
        abline(b=notch.large.D1_97.5, a=intersect.large.D1_97.5,col='grey')
        abline(b=notch.large.D1_2.5, a=intersect.large.D1_2.5,col='grey')

    plot(DF$fsc, DF$d2, xlim=c(0,65000), ylim=c(0,65000),main=paste(ins))
      points(df.c$fsc.med+comp.fsc,df.c$d2.med+comp.d2, col=2)
      points(ref.o$fsc.med, ref.o$d2.med, pch=16)
      points(ref.c$fsc.med+comp.fsc, ref.c$d2.med+comp.d2, col=2,pch=16)
      abline(b=notch.small.D2, a=intersect.small.D2,col=4)
      abline(b=notch.large.D2, a=intersect.large.D2,col=3)
        abline(b=notch.small.D2_97.5, a=intersect.small.D2_97.5,col='grey')
        abline(b=notch.small.D2_2.5, a=intersect.small.D2_2.5,col='grey')
        abline(b=notch.large.D2_97.5, a=intersect.large.D2_97.5,col='grey')
        abline(b=notch.large.D2_2.5, a=intersect.large.D2_2.5,col='grey')
}
dev.off()

write.csv(SLOPES, "/1.bead_calibration/seaflow_filter_MERGED_slopes.csv", quote=F, row.names=F)
