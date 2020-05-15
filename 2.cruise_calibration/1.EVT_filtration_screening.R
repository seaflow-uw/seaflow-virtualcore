library(popcycle)
library(dplyr)
library(splancs)
library(scales)

#########################################
## 0. DOWNLOAD THE RAW DATA using DAT ###
#########################################
# dat://fa2bb5981465ff583e091071f2f2c5e5f1b118c5219971a2dde107e4dec25e63

#Path to the raw data (DAT)
path.to.data <- "~/Documents/DATA/Codes/seaflow-virtualcore/seaflow-virtualcore-data/"
path.to.data <- "~/Documents/DATA/Codes/seaflow-virtualcore/"




####################################
### 1. SET FILTRATION PARAMETERS ###
####################################
beads.coord <- read.csv("https://raw.githubusercontent.com/armbrustlab/seaflow-filter/master/ALL-filterparams.csv")


#Path to the Git repository
setwd("~/Documents/DATA/Codes/seaflow-virtualcore/2.cruise_calibration/")


gate <- F
allcruises <- c("SCOPE_6", "SCOPE_2", "DeepDOM", "MBARI_1","SCOPE_16","Thompson_9", "MGL1704", "HOT")


for(cruise in allcruises){

  #cruise <- allcruises[8]

  print(cruise)
  ### Get EVT list
  if(cruise == "SCOPE_6" | cruise == "SCOPE_16" | cruise == "SCOPE_2"){
    list <- list.files(paste0(path.to.data,cruise,"data"), "00-00$",full.names=T, recursive=T)}
  if(cruise == "DeepDOM"| cruise == "MBARI_1" | cruise == "Thompson_9"){
    list <- list.files(paste0(path.to.data,cruise,"data"), ".evt$",full.names=T, recursive=T)}
  if(cruise == "Thompson_9"){
    list <- list.files(paste0(path.to.data,cruise,"data"), ".evt",full.names=T, recursive=T)}
  if(cruise == "SCOPE_2" | cruise == "MGL1704" | cruise == "HOT"){
    list <- list.files(paste0(path.to.data,cruise,"data"), "00-00.gz",full.names=T, recursive=T)}

  ### get Filtration parameters 
  # SLOPES
  slopes <- read.csv("seaflow_filter_slopes.csv")
  # BEADS
  if(cruise != "HOT"){
    beads.fsc <- as.numeric(beads.coord[which(beads.coord$cruise == cruise & beads.coord$quantile == 50),"beads.fsc.small"])
    beads.d1 <- as.numeric(beads.coord[which(beads.coord$cruise == cruise & beads.coord$quantile == 50),"beads.D1"])
    beads.d2 <- as.numeric(beads.coord[which(beads.coord$cruise == cruise & beads.coord$quantile == 50),"beads.D2"])
    # get instrument serial number
    inst <- as.numeric(beads.coord[which(beads.coord$cruise == cruise & beads.coord$quantile == 50),"instrument"])
  }

  # Alignment width
  width <- 5000

  ALL <- NULL
  for (file in list){
  #file <- list[2]
  print(file)

  evt <- readSeaflow(file,transform=F)

  # 1. Filtering noise
  evt. <- evt[evt$fsc_small > 1 | evt$D1 > 1 | evt$D2 > 1, ]

  # 2. Filtering out particles with saturated D1 and D2 signals
  evt.. <- evt.[evt.$D1 < max(evt.$D1) | evt.$D2 < max(evt.$D2), ]

  # 3. Fltering aligned particles (D1 = D2)
  aligned <- subset(evt.., D2 < D1 + width & D1 < D2 + width)

    # if(nrow(aligned)> 100000){aligned. <- sample_n(aligned, 100000)
    #   } else aligned. <- aligned
    # if(nrow(evt.)> 100000){df <- sample_n(evt., 100000)
    #   } else df <- evt.
    # par(mfrow=c(2,2))
    # plot_cyt(df, "D1", "D2")
    # plot_cyt(aligned., "D1", "D2")
    # plot_cyt(df, "fsc_small", "D1")
    # plot_cyt(df, "fsc_small", "pe")
    
    
    if(cruise == "HOT"){
      c <- basename(dirname(dirname(file)))
        beads.fsc <- as.numeric(beads.coord[which(beads.coord$cruise == c & beads.coord$quantile == 50),"beads.fsc.small"])
        beads.d1 <- as.numeric(beads.coord[which(beads.coord$cruise == c & beads.coord$quantile == 50),"beads.D1"])
        beads.d2 <- as.numeric(beads.coord[which(beads.coord$cruise == c & beads.coord$quantile == 50),"beads.D2"])
        inst <- as.numeric(beads.coord[which(beads.coord$cruise == c & beads.coord$quantile == 50),"instrument"])

      }

  # 4. Filtering focused particles (D/FSC)  
  screening <- c(0, 2500, 5000, 7500, 10000, 1, 2, 3, 4)
  for(corr in screening){

    print(paste("corr=", corr))
    notch.small.D1 <- slopes[slopes$ins== inst,'notch.small.D1']
    notch.small.D2 <- slopes[slopes$ins== inst,'notch.small.D2']
    notch.large.D1 <- slopes[slopes$ins== inst,'notch.large.D1']
    notch.large.D2 <- slopes[slopes$ins== inst,'notch.large.D2']
    offset.small.D1 <- round(beads.d1 - notch.small.D1 * beads.fsc) + corr
    offset.small.D2 <- round(beads.d2 - notch.small.D2 * beads.fsc) + corr
    offset.large.D1 <- round(beads.d1 - notch.large.D1 * beads.fsc) + corr
    offset.large.D2 <- round(beads.d2 - notch.large.D2 * beads.fsc) + corr
    
    if(corr == 1){
      offset.small.D1 <- 0
      offset.small.D2 <- 0
    }
    
    if(corr == 2){
      offset.small.D1 <- 0
      offset.small.D2 <- 0
      notch.small.D1 <- beads.d1/beads.fsc
      notch.small.D2 <- beads.d2/beads.fsc
    }
    
    if(corr == 3){
      offset.small.D1 <- 0
      offset.small.D2 <- 0
      notch.small.D1 <- (beads.d1+width/2)/beads.fsc
      notch.small.D2 <- (beads.d1+width/2)/beads.fsc
      offset.large.D1 <- round(beads.d1 - notch.large.D1 * beads.fsc) + width/2
      offset.large.D2 <- round(beads.d2 - notch.large.D2 * beads.fsc) + width/2
    }
    
    if(corr == 4){
      offset.small.D1 <- 0
      offset.small.D2 <- 0
      notch.small.D1 <- (beads.d1-width/2)/beads.fsc
      notch.small.D2 <- (beads.d2-width/2)/beads.fsc
      offset.large.D1 <- round(beads.d1 - notch.large.D1 * beads.fsc) - width/2
      offset.large.D2 <- round(beads.d2 - notch.large.D2 * beads.fsc) - width/2
    }
    
    # subset focused particles
    opp <- subset(aligned, D1 <= fsc_small*notch.small.D1 + offset.small.D1 & D2 <= fsc_small*notch.small.D2 + offset.small.D2 |
        D1  <= fsc_small*notch.large.D1 + offset.large.D1 & D2 <= fsc_small*notch.large.D2 + offset.large.D2)


      if(nrow(opp) < 10) next

      opp$pop <- 0
      opp. <- subset(opp, fsc_small > 2000)

        if(gate) {plot_cyt(opp , "fsc_small", "pe")
                        poly.beads <- getpoly(quiet=TRUE)
                        write.csv(poly.beads,paste0(cruise,"data/", basename(file),"-polybeads.csv"),quote=F, row.names=F)
                      }
        poly.beads <-  read.csv(paste0(cruise,"data/", basename(file),"-polybeads.csv"))
        b <- subset(opp,inout(opp[,c("fsc_small", "pe")],poly=poly.beads, bound=TRUE, quiet=TRUE))
        opp[row.names(b),'pop'] <- "beads"

        if(gate) {poly.syn <- getpoly(quiet=TRUE)
                        write.csv(poly.syn,paste0(cruise,"data/", basename(file),"-polysyn.csv"),quote=F, row.names=F)
                      }
        poly.syn <-  read.csv(paste0(cruise,"data/", basename(file),"-polysyn.csv"))
        s <- subset(opp,inout(opp[,c("fsc_small", "pe")],poly=poly.syn, bound=TRUE, quiet=TRUE))
        opp[row.names(s),'pop'] <- "synecho"

      x <- subset(opp, pop==0)
        if(gate){plot_cyt(opp , "fsc_small", "chl_small")
                        poly.pro <- getpoly(quiet=TRUE)
                        write.csv(poly.pro,paste0(cruise,"data/", basename(file),"-polypro.csv"),quote=F, row.names=F)
                      }
        poly.pro <-  read.csv(paste0(cruise,"data/", basename(file),"-polypro.csv"))
        p <- subset(x,inout(x[,c("fsc_small", "chl_small")],poly=poly.pro, bound=TRUE, quiet=TRUE))
        opp[row.names(p),'pop'] <- "prochloro"

      x <- subset(opp, pop==0)
      if(gate) { poly.pico <- getpoly(quiet=TRUE)
                      write.csv(poly.pico,paste0(cruise,"data/", basename(file),"-polypico.csv"),quote=F, row.names=F )
                    }
      poly.pico <-  read.csv(paste0(cruise,"data/", basename(file),"-polypico.csv"))
      l <- subset(x,inout(x[,c("fsc_small", "chl_small")],poly=poly.pico, bound=TRUE, quiet=TRUE))
      opp[row.names(l),'pop'] <- "pico"



    if(corr == 0 | corr == 1 | corr == 2 | corr == 3){
          png(paste0(cruise,"data/", basename(file),"-offset",corr,".png"),width=9, height=12, unit='in', res=100)

          if(nrow(opp)> 100000){opp. <- sample_n(opp, 100000)
             } else opp. <- opp
          par(mfrow=c(2,2))
          plot_cyt(aligned. ,  "fsc_small", "D1"); abline(b=notch.small.D1, a=offset.small.D1,col=2); abline(b=notch.large.D1, a=offset.large.D1,col=3); points(beads.fsc, beads.d1,pch=16, col=3)
          plot_cyt(aligned. , "fsc_small", "D2"); abline(b=notch.small.D2, a=offset.small.D2,col=2); abline(b=notch.large.D2, a=offset.large.D2,col=3); points(beads.fsc, beads.d2,pch=16, col=3)
          plot(opp.[, c("fsc_small", "pe")], xlim=c(0,66000), ylim=c(0,66000), cex = 1, pch = 16, col = alpha(as.numeric(as.factor(opp.$pop)), 0.4))
          plot(opp.[, c("fsc_small", "chl_small")], xlim=c(0,66000), ylim=c(0,66000), cex = 1, pch = 16, col = alpha(as.numeric(as.factor(opp.$pop)), 0.4))
          dev.off()
      }

      n.opp <- nrow(opp)
      n.opp. <- nrow(opp.)
      n.evt <- nrow(evt)
      n.evt. <- nrow(evt.)
      beads <- nrow(b)
      syn <- nrow(s)
      pro <- nrow(p)
      pico <- nrow(l)
      filename <- paste(basename(dirname(file)),basename(file),sep="/")
      all <- data.frame(cbind(file=filename, n.opp, n.opp., n.evt, n.evt., beads, pro, syn, pico, width, corr,
                                    offset.small.D1, offset.small.D2, offset.large.D1, offset.large.D2,
                                    notch.small.D1, notch.small.D2, notch.large.D1, notch.large.D2))

      ALL <- rbind(ALL, all)

    }
  }

write.csv(ALL, paste0(cruise,"data/seaflow-summary.csv"), quote=F, row.names=F)


}















##################
### 2. SUMMARY ###
##################
library(lattice)
library(lmodel2)
library(popcycle)

setwd("~/Documents/DATA/Codes/seaflow-virtualcore/2.cruise_calibration/")

allcruises <- c("SCOPE_6", "SCOPE_2", "DeepDOM", "MBARI_1","SCOPE_16","Thompson_9", "MGL1704", "HOT")
cruise <- allcruises[8]

DF <- NULL
for(cruise in allcruises){

print(cruise)

### REFERENCE ABUNDANCE from INFLUX
if(cruise == 'DeepDOM'){influx <- read.csv(paste0(cruise,"data/surface_samples_metadata.csv"))
                          influx$time <- as.POSIXct(influx$Date.Time, format="%m/%d/%y %H:%M", tz="GMT")
                        }
if(cruise == 'SCOPE_16'){influx <- read.csv(paste0(cruise,"data/all_station_curated.csv"))
                            influx$time_HST <- as.POSIXct(influx$time_HST, format="%m/%d/%y %H:%M", tz="HST")
                            influx <- influx[order(influx$time_HST), ]
                            influx <- subset(influx, Depth_m < 16)
                            pro <- subset(influx, population == 'prochloro')
                            syn <- subset(influx, population == 'synecho')
                            picoeuk <- subset(influx, population == 'picoeuk')
                          influx <- data.frame(cbind(time=pro$time_HST, pro=pro$abundance, syn=syn$abundance, pico=picoeuk$abundance))
                            influx$time <- as.POSIXct(influx$time, origin='1970-01-01',tz="GMT")
                              }
if(cruise == 'MBARI_1'){  sfl <- read.delim(paste0(cruise,"data/cruise.sfl"))
                          sfl$date <- as.POSIXct(sfl$DATE, format = "%FT%T", tz = "GMT")
                          inst <- 989
                          sfl$flow_rate2 <- flowrate(sfl$STREAM.PRESSURE, inst =inst)[,"flow_rate"]
                          influx <- read.csv(paste0(cruise,"data/mbari_summary.csv"))
                          id <- match(influx[,'file'],sfl[,"FILE"])
                          influx$time <- sfl[id, 'date']
                          influx$pro <- influx$pro / 1000
                          influx$syn <- influx$syn / 1000
                          influx$pico <- influx$picoeuks / 1000
                          }
if(cruise == 'SCOPE_6'){influx <- read.csv(paste0(cruise,"data/influx.csv"))
                          influx$time <- as.POSIXct(influx$time, tz="GMT")
                          pro <- subset(influx, pop == 'prochloro')
                          syn <- subset(influx, pop == 'synecho')
                          picoeuk <- subset(influx, pop == 'picoeuk')
                          crocco <- subset(influx, pop == 'crocco')
                        influx <- data.frame(cbind(time=pro$time, pro=pro$abundance, syn=syn$abundance, pico=picoeuk$abundance+crocco$abundance))
                          influx$time <- as.POSIXct(influx$time, origin='1970-01-01',tz="GMT")
                        }

if(cruise == 'Thompson_9'){pre.influx <- read.csv(paste0(cruise,"data/merged_summary_sds.csv"))
                            pre.influx$time <- as.POSIXct(pre.influx$local.time, format="%m/%d/%y %H:%M", tz="HST")
                            influx <- aggregate(pre.influx, by=list(pre.influx$fcs.file, pre.influx$i), FUN=mean)
                              pro <- subset(influx, Group.2 == 'prochloro')
                              syn <- subset(influx, Group.2 == 'synecho')
                            influx <- data.frame(cbind(time=pro$time, pro=pro$conc, syn=syn$conc, pico=NA))
                              influx$time <- as.POSIXct(influx$time, origin='1970-01-01',tz="GMT")
                        }


if(cruise == 'SCOPE_2'){influx <- read.csv(paste0(cruise,"data/abundance_table.csv"))
                          influx$time <- as.POSIXct(influx$time, format="%m/%d/%y %H:%M", tz="HST")
                          influx$pro <- influx[,"pro_infl"]
                          influx$syn <- influx[,"syn_infl"]
                          influx$pico <- influx[,"pico_infl"]
                        }

if(cruise == 'MGL1704'){influx <- read.csv(paste0(cruise,"data/all_station_curated.csv"))
                          influx$time <- as.POSIXct(influx$time_GMT, format="%m/%d/%y %H:%M", tz="GMT")
                          influx <- subset(influx, depth_m < 20)
                          pro <- subset(influx, population == 'prochloro')
                          syn <- subset(influx, population == 'synecho')
                          picoeuk <- subset(influx, population == 'picoeuk')
                          influx <- data.frame(cbind(time=pro$time, pro=pro$abundance, syn=syn$abundance, pico=picoeuk$abundance))
                          influx$time <- as.POSIXct(influx$time, origin='1970-01-01',tz="GMT")
                                        }




if(cruise == 'HOT'){influx <- read.csv(paste0(cruise,"data/HOT_influx_SeaFlow_results.csv"))
                          influx$time <- as.POSIXct(influx$HOT_datetime, format="%m/%d/%y %H:%M", tz="GMT")
                          influx <- data.frame(cbind(time=influx$time, pro=influx$HOT_abund_pro, syn=influx$HOT_abund_syn, pico=influx$HOT_abund_peuk))
                          influx$time <- as.POSIXct(influx$time, origin='1970-01-01',tz="GMT")
                                        }

influx <- influx[order(influx$time), ]

### SFL

if(cruise != 'MBARI_1' & cruise !="HOT"){
      inst <- as.numeric(beads.coord[which(beads.coord$cruise == cruise & beads.coord$quantile == 50),"instrument"])
      sfl <- read.delim(paste0("https://raw.githubusercontent.com/armbrustlab/seaflow-sfl/master/curated/", cruise,"_",inst, ".sfl"))
      sfl$date <- as.POSIXct(sfl$DATE, format="%FT%T", tz='GMT')
      sfl$flow_rate2 <- flowrate(sfl$STREAM.PRESSURE, inst =inst)[,"flow_rate"]
      }
      sfl$date <- as.POSIXct(sfl$DATE, format="%FT%T", tz='GMT')

if(cruise == "HOT"){
            list <- list.files(paste0(path.to.data,cruise,"data"), "00-00.gz",full.names=T, recursive=T)
              c <- basename(dirname(dirname(list)))
              sfl <- NULL
                for(i in c){
                  inst <- as.numeric(beads.coord[which(beads.coord$cruise == i & beads.coord$quantile == 50),"instrument"])
                  s <- read.delim(paste0("https://raw.githubusercontent.com/armbrustlab/seaflow-sfl/master/curated/", i,"_",inst, ".sfl"))
                  s$date <- as.POSIXct(s$DATE, format="%FT%T", tz='GMT')
                  s$flow_rate2 <- flowrate(s$STREAM.PRESSURE, inst =inst)[,"flow_rate"]
                  sfl <- rbind(sfl, s)
                  }
            }


### SEAFLOW DATA
ALL <- read.csv(paste0(cruise,"data/seaflow-summary.csv"))
if(cruise == "SCOPE_16") ALL <- ALL[!(ALL$file =="2016-04-26T15-07-38-00-00"),]
if(cruise == "DeepDOM") ALL <- ALL[!(ALL$file =="2013_094/321.evt" | ALL$file=="2013_124/409.evt"),]
if(cruise == "Thompson_9" | cruise == "SCOPE_2" | cruise == 'MGL1704' | cruise == "HOT") ALL$file <- sub(".gz","", ALL$file)

  ALL$cruise <- cruise
  id2 <- match(ALL[,'file'],as.character(sfl[,"FILE"]))

  ALL$time <- sfl[id2, 'date'] # add time
    ALL <- ALL[order(ALL$time), ]

  ALL$fr <- sfl[id2, 'flow_rate2'] # add flow rate

  # add abundance from INFLUX
  id <- findInterval(ALL$time,influx$time)
  if(cruise =="DeepDOM"| cruise == "SCOPE_16" | cruise == "Thompson_9"| cruise == "SCOPE_2" | cruise == "HOT") id <- id +1
  for (phyto in c('pro','syn','pico'))   ALL[,paste0(phyto,'.influx')] <- influx[id,phyto]



  ############### TEST different VC method ###################################
  correction <- 1
  if(cruise == "Thompson_9") correction <- 0.14 # correction factor due to troncated D1 / D2 (based on Seaflow_uncorrected abundance / influx)

  ALL$vc1 <- correction * 1000* 3*ALL[,'fr']*ALL[,'n.opp']/ALL[,'n.evt'] # calculate virtual core v1
  ALL$vc2 <- correction * 1000* 3*ALL[,'fr']*ALL[,'n.opp']/ALL[,'n.evt.'] # calculate virtual core v2
  ALL$vc3 <- correction * 1000* 3*ALL[,'fr']*ALL[,'n.opp.']/ALL[,'n.evt'] # calculate virtual core v3
  ALL$vc4 <- correction * 1000* 3*ALL[,'fr']*ALL[,'n.opp.']/ALL[,'n.evt.'] # calculate virtual core v4



  # calculate abundance from Seaflow
  for (phyto in c('pro','syn','pico')){
      ALL[,paste0(phyto,'.seaflow.each1')] <- ALL[,phyto] / ALL$vc1
      ALL[,paste0(phyto,'.seaflow.median1')] <- ALL[,phyto] / median(ALL$vc1)
      ALL[,paste0(phyto,'.seaflow.each2')] <- ALL[,phyto] / ALL$vc2
      ALL[,paste0(phyto,'.seaflow.median2')] <- ALL[,phyto] / median(ALL$vc2)
      ALL[,paste0(phyto,'.seaflow.each3')] <- ALL[,phyto] / ALL$vc3
      ALL[,paste0(phyto,'.seaflow.median3')] <- ALL[,phyto] / median(ALL$vc3)
      ALL[,paste0(phyto,'.seaflow.each4')] <- ALL[,phyto] / ALL$vc4
      ALL[,paste0(phyto,'.seaflow.median4')] <- ALL[,phyto] / median(ALL$vc4)
      }

      if(cruise == "Thompson_9") for(s in 1:4) ALL[,paste0('pro.seaflow.each',s)] <- ALL[,paste0('pro.seaflow.each',s)] + 150 # correction factor due to troncated D1 / D2 (based on Seaflow_corrected abundance - influx)

          # par(mfrow=c(3,1),cex=1.2)
          # for (phyto in c('pro','syn','pico')){
          #         plot(ALL[,'time'],ALL[,paste0(phyto,'.influx')],ylim=c(0,2*max(ALL[,paste0(phyto,'.influx')])),type='o', ylab="abundance", xlab=NA, main=paste(phyto))
          #         for(s in 1:4){
          #           points(ALL[,'time'],ALL[,paste0(phyto,'.seaflow.each',s)],col=s)
          #           points(ALL[,'time'],ALL[,paste0(phyto,'.seaflow.median',s)],col=s, pch=2)
          #           }
          #         }

    # calculate error between seaflow-based and influx-based abundance
    for (phyto in c('pro','syn','pico')){
      for(s in 1:4){
        ALL[,paste0(phyto,".diff.each",s)] <- abs(ALL[,paste0(phyto,'.seaflow.each',s)] - ALL[,paste0(phyto,'.influx')]) / ALL[,paste0(phyto,'.influx')]
        ALL[,paste0(phyto,".diff.median",s)] <- abs(ALL[,paste0(phyto,'.seaflow.median',s)] -ALL[,paste0(phyto,'.influx')]) / ALL[,paste0(phyto,'.influx')]
            }
      }

          # par(mfrow=c(3,1),cex=1.2)
          # for (phyto in c('pro','syn','pico')){
          #         plot(ALL[,'time'],rep(0, nrow(ALL)),ylim=c(0,1),type='l', ylab="Diff", xlab=NA, main=paste(phyto))
          #         for(s in 1:4){
          #           points(ALL[,'time'],ALL[,paste0(phyto,'.diff.each',s)],col=s)
          #           points(ALL[,'time'],ALL[,paste0(phyto,'.diff.median',s)],col=s, pch=2)
          #           }
          #         }


    DF <- rbind(DF,ALL)
  }


write.csv(DF,paste0("SeaflowInflux_comparison.csv"), quote=F, row.names=F)















###########################################
### 3. file-based VC vs cruise-based VC ###
###########################################
library(tidyverse)

setwd("~/Documents/DATA/Codes/seaflow-virtualcore/2.cruise_calibration/")
DF <- read_csv("SeaflowInflux_comparison.csv")
DF$corr <- factor(as.character(DF$corr), levels=c("0","2500","5000","7500","10000","1","2","3","4"))
  
# split by method
pos <- "stack"
pos <- "identity"

DF %>%
  filter(cruise !="Thompson_9") %>%
  ggplot(aes((pro.influx - pro.seaflow.each4)/pro.influx)) +
  geom_histogram(position=pos, alpha=0.5, fill="red3") +
  geom_histogram(aes((pro.influx - pro.seaflow.median4)/pro.influx), position=pos, alpha=0.5) +
  geom_vline(xintercept=0, lty=2) +
  scale_x_continuous(limits=c(-5,5)) +
  facet_wrap(~ corr) 

DF %>%
  filter(cruise !="Thompson_9") %>%
  ggplot(aes((syn.influx - syn.seaflow.each4)/syn.influx)) +
  geom_histogram(position=pos, alpha=0.5, fill='red3') +
  geom_histogram(aes((syn.influx - syn.seaflow.median4)/syn.influx), position=pos, alpha=0.5) +
  geom_vline(xintercept=0, lty=2) +
  scale_x_continuous(limits=c(-5,5)) +
  facet_wrap(~ corr)

DF %>%
  filter(cruise !="Thompson_9") %>%
  ggplot(aes((pico.influx - pico.seaflow.each4)/pico.influx)) +
  geom_histogram(position=pos, alpha=0.5, fill='red3') +
  geom_histogram(aes((pico.influx - pico.seaflow.median4)/pico.influx), position=pos, alpha=0.5) +
  geom_vline(xintercept=0, lty=2) +
  scale_x_continuous(limits=c(-5,5)) +
  facet_wrap(~ corr)

##############################
### 5. BEST VC calculation ###
##############################
# n.opp = total opp n.opp. = fsc_small > 2000; n.evt = total evt; n.evt. = D1 > 1 | D2 > 1 | fsc > 1
# n.opp/n.evt' #  v1   raw count 
# n.opp/n.evt. # v2  wo electrical noise
# n.opp./n.evt # v3  wo small particles
# n.opp./n.evt. # v4 wo small particles nor electrical noise

# split by method
w <- 0.1
DF %>%
  filter(cruise !="Thompson_9") %>%
  ggplot(aes((pro.influx - pro.seaflow.each1)/pro.influx)) +
  geom_histogram(binwidth=w, alpha=0.35, fill=1) +
  geom_histogram(aes((pro.influx - pro.seaflow.each2)/pro.influx), binwidth=w, alpha=0.35, fill=2) +
  geom_histogram(aes((pro.influx - pro.seaflow.each3)/pro.influx), binwidth=w, alpha=0.35, fill=3) +
  geom_histogram(aes((pro.influx - pro.seaflow.each4)/pro.influx), binwidth=w, alpha=0.35, fill=4) +
  geom_vline(xintercept=0, lty=2) +
  scale_x_continuous(limits=c(-5,5)) +
  facet_wrap(~ corr) 
print(paste0("best VC method: vc", 2))  
 
DF %>%
  filter(cruise !="Thompson_9") %>%
  ggplot(aes((syn.influx - syn.seaflow.each1)/syn.influx)) +
  geom_histogram(binwidth=w, alpha=0.35, fill=1) +
  geom_histogram(aes((syn.influx - syn.seaflow.each2)/syn.influx), binwidth=w, alpha=0.35, fill=2) +
  geom_histogram(aes((syn.influx - syn.seaflow.each3)/syn.influx), binwidth=w, alpha=0.35, fill=3) +
  geom_histogram(aes((syn.influx - syn.seaflow.each4)/syn.influx), binwidth=w, alpha=0.35, fill=4) +
  geom_vline(xintercept=0, lty=2) +
  scale_x_continuous(limits=c(-5,5)) +
  facet_wrap(~ corr) 
print(paste0("best VC method: vc", 2, " > ", 4))  

DF %>%
  filter(cruise !="Thompson_9") %>%
  ggplot(aes((pico.influx - pico.seaflow.median1)/pico.influx)) +
  geom_histogram(binwidth=w, alpha=0.35, fill=1) +
  geom_histogram(aes((pico.influx - pico.seaflow.median2)/pico.influx), binwidth=w, alpha=0.35, fill=2) +
  geom_histogram(aes((pico.influx - pico.seaflow.median3)/pico.influx), binwidth=w, alpha=0.35, fill=3) +
  geom_histogram(aes((pico.influx - pico.seaflow.median4)/pico.influx), binwidth=w, alpha=0.35, fill=4) +
  geom_vline(xintercept=0, lty=2) +
  scale_x_continuous(limits=c(-5,5)) +
  facet_wrap(~ corr) 
print(paste0("best VC method: vc", 2, " > ", 4))  


##############################
### 5. BEST SLOPE + OFFSET ###
##############################
library(viridis)
group.colors <- c(unknown='grey', beads='red3', prochloro=viridis(4)[1],synecho=viridis(4)[2],picoeuk=viridis(4)[3], croco=viridis(4)[4])
levels(df2$cruise) <- c("KN210-04","HOT","CN11","MGL1704" ,"KOK1606","KM1502","KM1513","TN271")

DF %>%
  filter(cruise !="Thompson_9") %>%
  ggplot(aes((pro.influx - pro.seaflow.each2)/pro.influx)) +
  geom_histogram(binwidth=w, position=pos, alpha=0.75, fill=group.colors[3]) +
  geom_histogram(aes((syn.influx - syn.seaflow.each2)/syn.influx), binwidth=w, position=pos, alpha=0.75, fill=group.colors[4]) +
  geom_histogram(aes((pico.influx - pico.seaflow.median2)/pico.influx), binwidth=w, position=pos, alpha=0.75, fill=group.colors[5]) +
  geom_vline(xintercept=0, lty=2) +
  scale_x_continuous(limits=c(-5,5)) +
  facet_wrap(~ corr) 
print(paste0("best SLOPE-OFFSET method: corr=",2))  





###################
### 7. PLOTTING ###
###################
library(tidyverse)
library(viridis)
library(ggpubr)

df1 <- subset(DF, corr== 2)

df2 <- aggregate(df1, by=list(df1$time), FUN=function(x) mean(x, na.rm=T))
  df2b <- aggregate(df1, by=list(df1$time), FUN=function(x) x[1])
  df2$cruise <- df2b$cruise

#remove SeaFlow pico data for Thompson_9
df2[which(df2$cruise == 'Thompson_9'), 'pico.seaflow.median4'] <- NA
  

df2 %>%
  ggplot(aes(pro.influx - pro.seaflow.each2, fill=cruise)) + 
  geom_histogram(cex=6, pch=21, alpha=0.5) +
  geom_vline(xintercept=0) +
  #scale_fill_viridis_d(option = "D") +
  #labs(x="time", y= "diff Influx - SeaFlow") +
  theme_bw()  

df2 %>%
  ggplot(aes(time, pro.seaflow.each2)) + 
  geom_point(cex=6, pch=21, alpha=0.5, fill=group.colors[3]) +
  geom_point(aes(time, syn.seaflow.each2), pch=21, alpha=0.5, cex=6, fill=group.colors[4]) +
  geom_point(aes(time, pico.seaflow.median2), pch=21, alpha=0.5, cex=6, fill=group.colors[5]) +
  geom_point(aes(time, pro.influx), pch=21, alpha=0.5, cex=2, fill=group.colors[3]) +
  geom_point(aes(time, syn.influx), pch=21, alpha=0.5, cex=2, fill=group.colors[4]) +
  geom_point(aes(time, pico.influx), pch=21, alpha=0.5, cex=2, fill=group.colors[5]) +
  facet_wrap(~ cruise, scales="free_x") +
  scale_y_continuous(trans= 'log10') +
  labs(x="", y=expression(paste("Abundance (cells µL"^{-1},")"))) +
  theme_bw()   

ggsave("SeaFlowInflux-CRUISEcomparison.png", width=12, height=9)
ggsave("SeaFlowInflux-CRUISEcomparison.pdf", width=12, height=9)





##### DIFFERENCE
s <- c(df2$pro.seaflow.each2, df2$syn.seaflow.each2, df2$pico.seaflow.median2)
i <- c(df2$pro.influx, df2$syn.influx, df2$pico.influx)
population <- rep(c("prochloro","synecho","picoeuk"), each=nrow(df2))
df <- tibble(s,i, population)

# calculate Percentage difference
# df$error <- 2*(df$i - df$s)/(df$i + df$s)

# calculate Log difference
df$error <- df$i/df$s

fold <- c(2,10)
print(paste(round(100*length(which(abs(df$error) < fold[1]))/nrow(df),2), "% of estimates have less than ",fold[1],"fold difference"))
print(paste(round(100*length(which(abs(df$error) > fold[2]))/nrow(df),2), "% of estimates have more than ",fold[2],"fold difference"))

r <- round(cor(df$i,df$s,use='pairwise.complete.obs',  method= 'pearson'),2)

p1 <- df %>%
    ggplot(aes(x=s, y=i)) +
    geom_point(aes(x=s, y=i, fill=population), pch=21, alpha=0.5, show.legend=T, cex=2) +
    geom_abline(intercept=0, col=1, lwd=0.5, lty=2) +
    geom_smooth(method='glm',col='red3',lwd=0.5) +
    scale_y_continuous(trans= 'log10') +
    scale_x_continuous(trans= 'log10') +
    scale_fill_manual(values=group.colors, guide=F) +
    labs(x=expression(paste("SeaFlow (cells µL"^{-1},")")), y=expression(paste("Influx (cells µL"^{-1},")"))) +
    annotate("text", x = 0.2, y=300, label = paste("r =",r)) +
    annotate("text", x = 0.2, y=200, label = paste("n =",nrow(df))) +
    theme_bw()

df$error <- log10(df$error)
p2 <- df %>%
    ggplot(aes(x=error)) + geom_histogram(aes(x=error, fill=population), binwidth=0.05, color='black', alpha=0.75, size=0.25) +
    scale_fill_manual(values=group.colors) +
    geom_vline(xintercept=log10(2), lwd=0.5, lty=2) +
    geom_vline(xintercept=-log10(2), lwd=0.5, lty=2) +
    labs(x=expression(paste("Log difference")), y=expression(paste("Count"))) +
    theme_bw()

ggarrange(p1, p2, labels = c("a", "b"), ncol = 2, nrow = 1)

#ggsave("SeaFlowInflux-correlation-v2.pdf", width=8, height=4)
