library(popcycle)
library(dplyr)
cols <- colorRampPalette(c("blue4", "royalblue4", "deepskyblue3", "seagreen3", "yellow", "orangered2", "darkred"))




######################################
## DOWNLOAD THE RAW DATA using DAT ###
######################################
dat://fa2bb5981465ff583e091071f2f2c5e5f1b118c5219971a2dde107e4dec25e63

#Path to the raw data (DAT)
path.to.data <- "~/Documents/DATA/Codes/seaflow-virtualcore/seaflow-virtualcore-data/"




#################################
### SET FILTRATION PARAMETERS ###
#################################
#Path to the Git repository
setwd("~/Documents/DATA/Codes/seaflow-virtualcore")

s <- 1

# SLOPES
slopes <- read.csv("1.bead_calibration/seaflow_filter_MERGED_slopes.csv")


gate <- FALSE

allcruises <- c("SCOPE6", "DeepDOM", "MBARI1","Gradient1")

cruise <- allcruises[3]

for(cruise in allcruises){

  print(cruise)
  ### Get EVT list
   if(cruise == "SCOPE6" | cruise == "Gradient1") list <- list.files(paste0("/2.cruise_calibration/",cruise,"data"), "00-00$",full.names=T, recursive=T)
   if(cruise == "DeepDOM"| cruise == "MBARI1") list <- list.files(paste0("/2.cruise_calibration/",cruise,"data"), ".evt$",full.names=T, recursive=T)

    if(cruise == "SCOPE6") inst <- 740
    if(cruise == "DeepDOM" | cruise == "MBARI1") inst <- 989
    if(cruise == "Gradient1") inst <- 751

    # BEADS
      if(cruise == "SCOPE6"){
        stat <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/stats.csv"))
        stat$time <- as.POSIXct(stat$time, format = "%FT%T", tz = "GMT")
        beads <- subset(stat, pop == 'beads')
        beads$fsc_small <- (log10(beads$fsc_small)/3.5)*2^16# LOG TRANSFORM DATA

        beads.fsc <- 58000
        plot(beads$time, beads$fsc_small,ylim=c(0,2^16))
        abline(h=beads.fsc, col=1)
        print(beads.fsc)
      }

      if(cruise == "MBARI1"){
        stat <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/stats.csv"))
        stat$time <- as.POSIXct(stat$time,  tz = "GMT")
        beads <- subset(stat, pop == 'beads')
        beads.fsc <- min(beads$fsc_small)
        plot(beads$time, beads$fsc_small,ylim=c(0,2^16))
        abline(h=beads.fsc, col=2)
        print(beads.fsc)
      }

      if(cruise == "DeepDOM" | cruise == "Gradient1"){
        stat <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/stats.csv"))
        stat$time <- as.POSIXct(stat$time, format = "%FT%T", tz = "GMT")
        beads <- subset(stat, pop == 'beads')
        beads$fsc_small <- (log10(beads$fsc_small)/3.5)*2^16# LOG TRANSFORM DATA

        beads.fsc <- median(beads$fsc_small)
        plot(beads$time, beads$fsc_small,ylim=c(0,2^16))
        abline(h=beads.fsc, col=2)
        print(beads.fsc)
      }

      notch.small.D1ref <- slopes[slopes$ins== inst,'notch.small.D1']
      notch.small.D2ref <- slopes[slopes$ins== inst,'notch.small.D2']

    beads <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/beadsD1D2.csv"))
      beads.d1 <- median(beads$beads.d1);print(beads.d1)
      beads.d2 <- median(beads$beads.d2);print(beads.d2)
        correction.D1 <- round(beads.d1 - notch.small.D1ref * beads.fsc)
        correction.D2 <- round(beads.d2 - notch.small.D2ref * beads.fsc)


  width <- 2500
  ALL <- NULL

    for (file in list){

      #file <- list[4]
      print(file)

      evt <- readSeaflow(file,transform=F)

      # Filtering noise
      evt. <- evt[evt$fsc_small > 1 | evt$D1 > 1 | evt$D2 > 1, ]

      # Correction for the difference of sensitivity between D1 and D2
      origin <- 0#median(evt.[evt.$D2>5000,'D2'])-median(evt.[evt.$D1>5000,'D1'])

      # Fltering aligned particles (D1 = D2), with Correction for the difference of sensitivity between D1 and D2
      aligned <- subset(evt., D2 < (D1+origin) + width & (D1+origin) < D2 + width )


        if(nrow(aligned)> 100000){aligned. <- sample_n(aligned, 100000)
          } else aligned. <- aligned


       for(i in c("mean", "2.5","97.5")){

              print(i)

                if(i == "mean"){
                  notch.small.D1 <- slopes[slopes$ins== inst,'notch.small.D1']
                  notch.small.D2 <- slopes[slopes$ins== inst,'notch.small.D2']
                  notch.large.D1 <- slopes[slopes$ins== inst,'notch.large.D1']
                  notch.large.D2 <- slopes[slopes$ins== inst,'notch.large.D2']
                  intersect.small.D1 <- slopes[slopes$ins== inst,'intersect.small.D1']
                  intersect.small.D2 <- slopes[slopes$ins== inst,'intersect.small.D2']
                  intersect.large.D1 <- slopes[slopes$ins== inst,'intersect.large.D1']
                  intersect.large.D2 <- slopes[slopes$ins== inst,'intersect.large.D2']
                     }
                if(i == "97.5"){
                  notch.small.D1 <- slopes[slopes$ins== inst,'notch.small.D1_97.5']
                  notch.small.D2 <- slopes[slopes$ins== inst,'notch.small.D2_97.5']
                  notch.large.D1 <- slopes[slopes$ins== inst,'notch.large.D1_97.5']
                  notch.large.D2 <- slopes[slopes$ins== inst,'notch.large.D2_97.5']
                  intersect.small.D1 <- slopes[slopes$ins== inst,'intersect.small.D1_97.5']
                  intersect.small.D2 <- slopes[slopes$ins== inst,'intersect.small.D2_97.5']
                  intersect.large.D1 <- slopes[slopes$ins== inst,'intersect.large.D1_97.5']
                  intersect.large.D2 <- slopes[slopes$ins== inst,'intersect.large.D2_97.5']
                  }
                if(i == "2.5"){
                  notch.small.D1 <- slopes[slopes$ins== inst,'notch.small.D1_2.5']
                  notch.small.D2 <- slopes[slopes$ins== inst,'notch.small.D2_2.5']
                  notch.large.D1 <- slopes[slopes$ins== inst,'notch.large.D1_2.5']
                  notch.large.D2 <- slopes[slopes$ins== inst,'notch.large.D2_2.5']
                  intersect.small.D1 <- slopes[slopes$ins== inst,'intersect.small.D1_2.5']
                  intersect.small.D2 <- slopes[slopes$ins== inst,'intersect.small.D2_2.5']
                  intersect.large.D1 <- slopes[slopes$ins== inst,'intersect.large.D1_2.5']
                  intersect.large.D2 <- slopes[slopes$ins== inst,'intersect.large.D2_2.5']
                }


        ### OFFSET for BEST 1.0
        if(s==1){
          corr <- 0
          offset.small.D1 <- round(beads.d1 - notch.small.D1 * beads.fsc - correction.D1) + corr
          offset.small.D2 <- round(beads.d2 - notch.small.D2 * beads.fsc - correction.D2) + corr
          offset.large.D1 <- round(beads.d1 - notch.large.D1 * beads.fsc - correction.D1) + corr
          offset.large.D2 <- round(beads.d2 - notch.large.D2 * beads.fsc - correction.D2) + corr
          }

        ### OFFSET for BEST 2.0
        if(s == 2){
          offset.small.D1 <- 0
          offset.small.D2 <- 0
          offset.large.D1 <- intersect.large.D1 - intersect.small.D1
          offset.large.D2 <- intersect.large.D2 - intersect.small.D2
          }


      opp <- subset(aligned, D1 <= fsc_small*notch.small.D1 + offset.small.D1 & D2 <= fsc_small*notch.small.D2 + offset.small.D2 |
                            D1 <= fsc_small*notch.large.D1 + offset.large.D1 & D2 <= fsc_small*notch.large.D2 + offset.large.D2)

          if(nrow(opp) < 10) next

      opp$pop <- 0

      print("Beads")
        if(gate) {plot.cytogram(opp , "fsc_small", "pe")
                        poly.beads <- getpoly(quiet=TRUE)
                        write.csv(poly.beads,paste0("/2.cruise_calibration/",cruise,"data/", basename(file),"-polybeads.csv"),quote=F, row.names=F)
                      }
        poly.beads <-  read.csv(paste0("/2.cruise_calibration/",cruise,"data/", basename(file),"-polybeads.csv"))
        b <- subset(opp,inout(opp[,c("fsc_small", "pe")],poly=poly.beads, bound=TRUE, quiet=TRUE))
        opp[row.names(b),'pop'] <- "beads"

      print("Syn")
        if(gate) {poly.syn <- getpoly(quiet=TRUE)
                        write.csv(poly.syn,paste0("/2.cruise_calibration/",cruise,"data/", basename(file),"-polysyn.csv"),quote=F, row.names=F)
                      }
        poly.syn <-  read.csv(paste0("/2.cruise_calibration/",cruise,"data/", basename(file),"-polysyn.csv"))
        a <- subset(opp,inout(opp[,c("fsc_small", "pe")],poly=poly.syn, bound=TRUE, quiet=TRUE))
        opp[row.names(a),'pop'] <- "synecho"

      print("Pro")
      x <- subset(opp, pop==0)
        if(gate){plot.cytogram(opp , "fsc_small", "chl_small")
                        poly.pro <- getpoly(quiet=TRUE)
                        write.csv(poly.pro,paste0("/2.cruise_calibration/",cruise,"data/", basename(file),"-polypro.csv"),quote=F, row.names=F)
                      }
        poly.pro <-  read.csv(paste0("/2.cruise_calibration/",cruise,"data/", basename(file),"-polypro.csv"))
        p <- subset(x,inout(x[,c("fsc_small", "chl_small")],poly=poly.pro, bound=TRUE, quiet=TRUE))
        opp[row.names(p),'pop'] <- "prochloro"

      print("Pico")
      x <- subset(opp, pop==0)
        if(gate) { poly.pico <- getpoly(quiet=TRUE)
                        write.csv(poly.pico,paste0("/2.cruise_calibration/",cruise,"data/", basename(file),"-polypico.csv"),quote=F, row.names=F )
                      }
        poly.pico <-  read.csv(paste0("/2.cruise_calibration/",cruise,"data/", basename(file),"-polypico.csv"))
        l <- subset(x,inout(x[,c("fsc_small", "chl_small")],poly=poly.pico, bound=TRUE, quiet=TRUE))
        opp[row.names(l),'pop'] <- "pico"



          if(i == "mean"){
          png(paste0("/2.cruise_calibration/",cruise,"data/", basename(file),"-corr",corr,".png"),width=9, height=12, unit='in', res=100)
          if(nrow(opp)> 100000){opp. <- sample_n(opp, 100000)
             } else opp. <- opp
          par(mfrow=c(2,2))
          plot.cytogram(aligned. ,  "fsc_small", "D1"); abline(b=notch.small.D1, a=offset.small.D1,col=2); abline(b=notch.large.D1, a=offset.large.D1,col=3)
          plot.cytogram(aligned. , "fsc_small", "D2"); abline(b=notch.small.D2, a=offset.small.D2,col=2); abline(b=notch.large.D2, a=offset.large.D2,col=3)
          plot.vct.cytogram(opp. , "fsc_small", "pe")
          plot.vct.cytogram(opp. , "fsc_small", "chl_small")

          dev.off()
          }

      n.opp <- nrow(opp)
      n.opp2 <- nrow(subset(opp, fsc_small > 1))
      n.evt <- nrow(evt)
      n.evt. <- nrow(evt.)
      beads <- nrow(b)
      syn <- nrow(a)
      pro <- nrow(p)
      pico <- nrow(l)
      filename <- basename(file)
      if(cruise == 'MBARI1' | cruise == 'DeepDOM') filename <- paste(basename(dirname(file)),basename(file),sep="/")

      all <- data.frame(cbind(file=filename, i, n.opp, n.opp2, n.evt, n.evt., beads, pro, syn, pico, origin, offset.small.D1, offset.large.D1,offset.small.D2, offset.large.D2, width))

      ALL <- rbind(ALL, all)

          }
      }

write.csv(ALL, paste0("/2.cruise_calibration/",cruise,"data/seaflow-summaryBEST",s,".csv"), quote=F, row.names=F)


}















###############
### SUMMARY ###
###############
library(lattice)
library(lmodel2)
library(plotrix)

s <- 2

#Path to the Git repository
setwd("~/Documents/DATA/Codes/seaflow-virtualcore")



allcruises <- c("SCOPE6", "DeepDOM", "MBARI1","Gradient1")
cruise <- allcruises[2]



summary <- NULL
for(cruise in allcruises){

print(cruise)

### SFL + INFLUX
if(cruise == 'DeepDOM'){sfl <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/sfl.csv"))
                        sfl$date <- as.POSIXct(sfl$date, format = "%FT%T", tz = "GMT")
                      influx <- read.delim(paste0("/2.cruise_calibration/",cruise,"data/surface_samples_metadata.tab"))
                        influx$time <- as.POSIXct(influx$Date.Time, format="%m/%d/%y %H:%M", tz="GMT")
                        }
if(cruise == 'Gradient1'){sfl <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/sfl.csv"))
                            sfl$date <- as.POSIXct(sfl$date, format = "%FT%T", tz = "GMT")
                            sfl$file <- basename(as.character(sfl$file))
                          influx <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/all_station_curated.csv"))
                            influx$time_HST <- as.POSIXct(influx$time_HST, format="%m/%d/%y %H:%M", tz="HST")
                            influx <- influx[order(influx$time_HST), ]
                            influx <- subset(influx, Depth_m < 16)
                            pro <- subset(influx, population == 'prochloro')
                            syn <- subset(influx, population == 'synecho')
                            picoeuk <- subset(influx, population == 'picoeuk')
                          influx <- data.frame(cbind(time=pro$time_HST, pro=pro$abundance, syn=syn$abundance, pico=picoeuk$abundance))
                            influx$time <- as.POSIXct(influx$time, origin='1970-01-01',tz="GMT")
                              }
if(cruise == 'MBARI1'){sfl <- read.delim(paste0("/2.cruise_calibration/",cruise,"data/cruise.sfl"))
                        sfl$date <- as.POSIXct(sfl$DATE, format = "%FT%T", tz = "GMT")
                        sfl$file <- sfl[,"FILE"]
                        sfl$flow_rate <- sfl[,"FLOW.RATE"]
                      influx <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/mbari_summary.csv"))
                        id <- match(influx[,'file'],sfl[,"FILE"])
                        influx$time <- sfl[id, 'date']
                        influx$pro <- influx$pro / 1000
                        influx$syn <- influx$syn / 1000
                        influx$pico <- influx$picoeuks / 1000
                          }
if(cruise == 'SCOPE6'){sfl <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/sfl.csv"))
                        sfl$date <- as.POSIXct(sfl$date, format = "%FT%T", tz = "GMT")
                        sfl$flow_rate <- median(sfl$flow_rate, na.rm=T)
                        sfl$file <- basename(as.character(sfl$file))
                      influx <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/influx.csv"))
                        influx$time <- as.POSIXct(influx$time, tz="GMT")
                        pro <- subset(influx, pop == 'prochloro')
                        syn <- subset(influx, pop == 'synecho')
                        picoeuk <- subset(influx, pop == 'picoeuk')
                        crocco <- subset(influx, pop == 'crocco')
                          influx <- data.frame(cbind(time=pro$time, pro=pro$abundance, syn=syn$abundance, pico=picoeuk$abundance+crocco$abundance))
                          influx$time <- as.POSIXct(influx$time, origin='1970-01-01',tz="GMT")
                        }
influx <- influx[order(influx$time), ]

### SEAFLOW
ALL <- read.csv(paste0("/2.cruise_calibration/",cruise,"data/seaflow-summaryBEST",s,".csv"))
  if(cruise == "Gradient1") ALL <- ALL[!(ALL$file =="2016-04-26T15-07-38-00-00"),]
  if(cruise == "DeepDOM") ALL <- ALL[!(ALL$file =="2013_094/321.evt" | ALL$file=="2013_124/409.evt"),]

  id2 <- match(ALL[,'file'],as.character(sfl[,"file"]))
  ALL$time <- sfl[id2, 'date']
  ALL$fr <- sfl[id2, 'flow_rate']
  ALL <- ALL[order(ALL$time), ]
  ALL$vc <- 3*ALL[,'fr']*ALL[,'n.opp2']/ALL[,'n.evt.']
  ALL[which(ALL$pro < 110),'pro'] <- 1

#     plot(ALL$time, ALL[,'n.opp'])
#     points(ALL$time, ALL[,'n.opp2'],col=2)


  ALL[,'pro.seaflow'] <- ALL[,'pro']/ALL[,'vc']
  ALL[,'pro.seaflow2'] <- ALL[,'pro']/median(ALL[,'vc'])
  ALL[,'syn.seaflow'] <- ALL[,'syn']/ALL[,'vc']
  ALL[,'syn.seaflow2'] <- ALL[,'syn']/median(ALL[,'vc'])
  ALL[,'pico.seaflow'] <- ALL[,'pico']/ALL[,'vc']
  ALL[,'pico.seaflow2'] <- ALL[,'pico']/median(ALL[,'vc'])


# plot(ALL$corr, ALL$vc,pch=2, log='y')

id <- findInterval(ALL$time,influx$time)
if(cruise =="DeepDOM"| cruise == "Gradient1") id <- id +1
ALL[,'pro.influx'] <- influx[id,'pro']
ALL[,'syn.influx'] <- influx[id,'syn']
ALL[,'pico.influx'] <- influx[id,'pico']

  # ### PLOTING ALL
  # par(mfrow=c(3,1),cex=1.2, pty='m')
  # plot(ALL$time, ALL[,'pro']/ALL[,'vc'], col=1, main="pro", ylab="Abundance", xlab=NA)
  #    points(ALL$time,ALL[,'pro.influx'],pch=1, cex=2, col='darkgrey',type='o')
  # plot(ALL$time, ALL[,'syn']/ALL[,'vc'], col=1, main="syn", ylab="Abundance", xlab=NA)
  #   points(ALL$time,ALL[,'syn.influx'],pch=1, cex=2, col='darkgrey',type='o')
  # plot(ALL$time, ALL[,'pico']/ALL[,'vc'], col=1, main="pico", ylab="Abundance", xlab=NA)
  #   points(ALL$time,ALL[,'pico.influx'],pch=1, cex=2, col='darkgrey',type='o')


png(paste0("/2.cruise_calibration/",cruise, "-SeaFlowInflux-compALL-BEST",s,".png"),width=12, height=12, unit='in', res=100)

    par(mfrow=c(3,2),cex=1.2, mar=c(4,4,2,1))

    df <- subset(ALL, i == "mean");    df1 <- subset(ALL, i == "2.5");    df2 <- subset(ALL, i == "97.5")

    par(pty='m')
    plot(df$time, df[,'pro']/df[,'vc'], main=paste("pro"), ylab="Abundance", xlab=NA,type='o')
      lines(df1$time, df1[,'pro']/df1[,'vc'],col='grey')
      lines(df2$time, df2[,'pro']/df2[,'vc'],col='grey')
    points(df$time, df[,'pro.influx'],pch=16, cex=2, col='darkgrey',type='o')

    par(pty='s')
    plot(df[,'pro']/df[,'vc'],df[,'pro.influx'],ylab="Influx", xlab="SeaFlow"); abline(b=1,a=0)


    par(pty='m')
    plot(df$time, df[,'syn']/median(df[,'vc']), main=paste("syn"), ylab="Abundance", xlab=NA,type='o')
      lines(df1$time, df1[,'syn']/median(df1[,'vc']),col='grey')
      lines(df2$time, df2[,'syn']/median(df2[,'vc']),col='grey')
    points(df$time, df[,'syn.influx'],pch=16, cex=2, col='darkgrey',type='o')

    par(pty='s')
    plot(df[,'syn']/median(df[,'vc']),df[,'syn.influx'],ylab="Influx", xlab="SeaFlow"); abline(b=1,a=0)


    par(pty='m')
    plot(df$time, df[,'pico']/median(df[,'vc']), main=paste("pico"), ylab="Abundance", xlab=NA,type='o')
      lines(df1$time, df1[,'pico']/median(df1[,'vc']),col='grey')
      lines(df2$time, df2[,'pico']/median(df2[,'vc']),col='grey')
    points(df$time, df[,'pico.influx'],pch=16, cex=2, col='darkgrey',type='o')

    par(pty='s')
    plot(df[,'pico']/median(df[,'vc']),df[,'pico.influx'],ylab="Influx", xlab="SeaFlow"); abline(b=1,a=0)

  dev.off()


        DF <- data.frame(cbind(cruise, df ))
  summary <- rbind(summary, DF)



  # results <- results[order(results$offset),]
  # results <- results[order(results$vc),]
  # par(mfrow=c(1,1))
  # plot(results$offset, results$diff.pro, log='y', ylim=c(min(results$diff.pico),max(results$diff.pro)),pch=as.numeric(results$vc),type='o')
  # points(results$offset, results$diff.syn, col=2, pch=as.numeric(results$vc),type='o')
  # points(results$offset, results$diff.pico, col=3,pch=as.numeric(results$vc),type='o')




}

write.csv(summary,paste0("/2.cruise_calibration/SeaFlow-Influx_compALL_BEST",s,".csv"), quote=F, row.names=F)






############
### BEST ###
############
library(lmodel2)
#Path to the Git repository
setwd("~/Documents/DATA/Codes/seaflow-virtualcore")
summary <- read.csv(paste0("/2.cruise_calibration/SeaFlow-Influx_compALL_BEST",s,".csv"))

png(paste0("/2.cruise_calibration/SeaFlowInflux-compALL_BEST",s,".png"),width=6, height=12, unit='in', res=100)
  par(mfrow=c(3,1),cex=1.2, mar=c(4,4,2,1),pty='s')
  plot(summary[,'pro.seaflow'],summary[,'pro.influx'],asp=1,ylab="Influx", xlab="SeaFlow",col=summary[,'cruise'], pch=16, ylim=c(50,500), xlim=c(50,500) ,log='xy'); abline(b=1,a=0,lwd=2,col='grey')
  plot(summary[,'syn.seaflow2'],summary[,'syn.influx'],asp=1,ylab="Influx", xlab="SeaFlow",col=summary[,'cruise'], pch=16, ylim=c(0.5,200), xlim=c(0.5,200), log='xy'); abline(b=1,a=0,lwd=2,col='grey')
  plot(summary[,'pico.seaflow2'],summary[,'pico.influx'],asp=1,ylab="Influx", xlab="SeaFlow",col=summary[,'cruise'], pch=16, ylim=c(0.5,100), xlim=c(0.5,100), log='xy'); abline(b=1,a=0,lwd=2,col='grey')
dev.off()
