library(popcycle)
library(dplyr)
cols <- colorRampPalette(c("blue4", "royalblue4", "deepskyblue3", "seagreen3", "yellow", "orangered2", "darkred"))


#########################################
## 0. DOWNLOAD THE RAW DATA using DAT ###
#########################################
# dat://fa2bb5981465ff583e091071f2f2c5e5f1b118c5219971a2dde107e4dec25e63

#Path to the raw data (DAT)
path.to.data <- "~/Documents/DATA/Codes/seaflow-virtualcore/seaflow-virtualcore-data/"




####################################
### 1. SET FILTRATION PARAMETERS ###
####################################
beads.coord <- read.csv("https://raw.githubusercontent.com/armbrustlab/seaflow-filter/master/ALL-filterparams.csv")


#Path to the Git repository
setwd("~/Documents/DATA/Codes/seaflow-virtualcore/2.cruise_calibration/")


gate <- FALSE
allcruises <- c("SCOPE_6", "DeepDOM", "MBARI_1","SCOPE_16","Thompson_9")

for(cruise in allcruises){

  #cruise <- allcruises[5]

  print(cruise)
  ### Get EVT list
  if(cruise == "SCOPE_6" | cruise == "SCOPE_16") list <- list.files(paste0(path.to.data,cruise,"data"), "00-00$",full.names=T, recursive=T)
  if(cruise == "DeepDOM"| cruise == "MBARI_1" | cruise == "Thompson_9") list <- list.files(paste0(path.to.data,cruise,"data"), ".evt$",full.names=T, recursive=T)
  if(cruise == "Thompson_9") list <- list.files(paste0(path.to.data,cruise,"data"), ".evt",full.names=T, recursive=T)


  if(cruise == "SCOPE_6" | cruise == "Thompson_9") inst <- 740
  if(cruise == "DeepDOM" | cruise == "MBARI_1") inst <- 989
  if(cruise == "SCOPE_16") inst <- 751

  # SLOPES
  slopes <- read.csv("seaflow_filter_slopes.csv")

  # BEADS
  beads.fsc <- as.numeric(beads.coord[which(beads.coord$cruise == cruise & beads.coord$quantile == 50),"beads.fsc.small"])
  beads.d1 <- as.numeric(beads.coord[which(beads.coord$cruise == cruise & beads.coord$quantile == 50),"beads.D1"])
  beads.d2 <- as.numeric(beads.coord[which(beads.coord$cruise == cruise & beads.coord$quantile == 50),"beads.D2"])

width <- 3000
ALL <- NULL

for (file in list){

  #file <- list[4]
  print(file)

  evt <- readSeaflow(file,transform=F)

  # Filtering noise
  evt. <- evt[evt$fsc_small > 1 | evt$D1 > 1 | evt$D2 > 1, ]

  # Fltering aligned particles (D1 = D2), with Correction for the difference of sensitivity between D1 and D2
  aligned <- subset(evt., D2 < D1 + width & D1 < D2 + width)

    if(nrow(aligned)> 100000){aligned. <- sample_n(aligned, 100000)
      } else aligned. <- aligned

    if(nrow(evt.)> 100000){df <- sample_n(evt., 100000)
      } else df <- evt.
    # par(mfrow=c(2,2))
    # plot.cytogram(df, "D1", "D2"); abline(b=notch, a=origin)
    # plot.cytogram(aligned., "D1", "D2"); abline(b=1, a=origin)
    # plot.cytogram(df, "fsc_small", "D1"); abline(b=notch.small.D1, a=offset.small.D1,col=2); abline(b=notch.large.D1, a=offset.large.D1,col=3); points(beads.fsc, beads.d1,pch=16, col=3)
    # plot.cytogram(df, "fsc_small", "D2"); abline(b=notch.small.D2, a=offset.small.D2,col=2); abline(b=notch.large.D2, a=offset.large.D2,col=3); points(beads.fsc, beads.d2,pch=16, col=3)

    screening <- c(-5000,-2500,0,2500,5000,7500,10000, 1, 2, 3, 4)

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
        notch.small.D1 <- (beads.d1+width)/beads.fsc
        notch.small.D2 <- (beads.d1+width)/beads.fsc
        offset.large.D1 <- round(beads.d1 - notch.large.D1 * beads.fsc) + width
        offset.large.D2 <- round(beads.d2 - notch.large.D2 * beads.fsc) + width
        }

    if(corr == 4){
        offset.small.D1 <- 0
        offset.small.D2 <- 0
        notch.small.D1 <- (beads.d1-width)/beads.fsc
        notch.small.D2 <- (beads.d2-width)/beads.fsc
        offset.large.D1 <- round(beads.d1 - notch.large.D1 * beads.fsc) - width
        offset.large.D2 <- round(beads.d2 - notch.large.D2 * beads.fsc) - width
        }

    opp <- subset(aligned, D1 <= fsc_small*notch.small.D1 + offset.small.D1 & D2 <= fsc_small*notch.small.D2 + offset.small.D2 |
        D1  <= fsc_small*notch.large.D1 + offset.large.D1 & D2 <= fsc_small*notch.large.D2 + offset.large.D2)


      if(nrow(opp) < 10) next

      opp$pop <- 0
      opp. <- subset(opp, fsc_small > 2000)

      print("Beads")
        if(gate) {plot.cytogram(opp , "fsc_small", "pe")
                        poly.beads <- getpoly(quiet=TRUE)
                        write.csv(poly.beads,paste0(cruise,"data/", basename(file),"-polybeads.csv"),quote=F, row.names=F)
                      }
        poly.beads <-  read.csv(paste0(cruise,"data/", basename(file),"-polybeads.csv"))
        b <- subset(opp,inout(opp[,c("fsc_small", "pe")],poly=poly.beads, bound=TRUE, quiet=TRUE))
        opp[row.names(b),'pop'] <- "beads"

      print("Syn")
        if(gate) {poly.syn <- getpoly(quiet=TRUE)
                        write.csv(poly.syn,paste0(cruise,"data/", basename(file),"-polysyn.csv"),quote=F, row.names=F)
                      }
        poly.syn <-  read.csv(paste0(cruise,"data/", basename(file),"-polysyn.csv"))
        s <- subset(opp,inout(opp[,c("fsc_small", "pe")],poly=poly.syn, bound=TRUE, quiet=TRUE))
        opp[row.names(s),'pop'] <- "synecho"

      print("Pro")
      x <- subset(opp, pop==0)
        if(gate){plot.cytogram(opp , "fsc_small", "chl_small")
                        poly.pro <- getpoly(quiet=TRUE)
                        write.csv(poly.pro,paste0(cruise,"data/", basename(file),"-polypro.csv"),quote=F, row.names=F)
                      }
        poly.pro <-  read.csv(paste0(cruise,"data/", basename(file),"-polypro.csv"))
        p <- subset(x,inout(x[,c("fsc_small", "chl_small")],poly=poly.pro, bound=TRUE, quiet=TRUE))
        opp[row.names(p),'pop'] <- "prochloro"

      print("Pico")
      x <- subset(opp, pop==0)
      if(gate) { poly.pico <- getpoly(quiet=TRUE)
                      write.csv(poly.pico,paste0(cruise,"data/", basename(file),"-polypico.csv"),quote=F, row.names=F )
                    }
      poly.pico <-  read.csv(paste0(cruise,"data/", basename(file),"-polypico.csv"))
      l <- subset(x,inout(x[,c("fsc_small", "chl_small")],poly=poly.pico, bound=TRUE, quiet=TRUE))
      opp[row.names(l),'pop'] <- "pico"



    if(corr == 0 | corr == 2 | corr == 3 | corr == 4){
          png(paste0(cruise,"data/", basename(file),"-offset",corr,".png"),width=9, height=12, unit='in', res=100)

          if(nrow(opp)> 100000){opp. <- sample_n(opp, 100000)
             } else opp. <- opp
          par(mfrow=c(2,2))
          plot.cytogram(aligned. ,  "fsc_small", "D1"); abline(b=notch.small.D1, a=offset.small.D1,col=2); abline(b=notch.large.D1, a=offset.large.D1,col=3); points(beads.fsc, beads.d1,pch=16, col=3)
          plot.cytogram(aligned. , "fsc_small", "D2"); abline(b=notch.small.D2, a=offset.small.D2,col=2); abline(b=notch.large.D2, a=offset.large.D2,col=3); points(beads.fsc, beads.d2,pch=16, col=3)
          plot.vct.cytogram(opp. , "fsc_small", "pe")
          plot.vct.cytogram(opp. , "fsc_small", "chl_small")

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

allcruises <- c("SCOPE_6", "DeepDOM", "MBARI_1","SCOPE_16","Thompson_9")
cruise <- allcruises[4]
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

influx <- influx[order(influx$time), ]


### SFL
if(cruise == "SCOPE_6" | cruise == "Thompson_9") inst <- 740
if(cruise == "DeepDOM" | cruise == "MBARI_1") inst <- 989
if(cruise == "SCOPE_16") inst <- 751

if(cruise != 'MBARI_1'){
      sfl <- read.delim(paste0("https://raw.githubusercontent.com/armbrustlab/seaflow-sfl/master/curated/", cruise,"_",inst, ".sfl"))
      sfl$date <- as.POSIXct(sfl$DATE, format="%FT%T", tz='GMT')
      }

### SEAFLOW DATA
ALL <- read.csv(paste0(cruise,"data/seaflow-summary.csv"))
if(cruise == "SCOPE_16") ALL <- ALL[!(ALL$file =="2016-04-26T15-07-38-00-00"),]
if(cruise == "DeepDOM") ALL <- ALL[!(ALL$file =="2013_094/321.evt" | ALL$file=="2013_124/409.evt"),]
if(cruise == "Thompson_9") ALL$file <- sub(".gz","", ALL$file)

  ALL$cruise <- cruise
  id2 <- match(ALL[,'file'],as.character(sfl[,"FILE"]))

  ALL$time <- sfl[id2, 'date'] # add time
    ALL <- ALL[order(ALL$time), ]

  sfl$flow_rate2 <- flowrate(sfl$STREAM.PRESSURE, inst =inst)[,"flow_rate"]
  ALL$fr <- sfl[id2, 'flow_rate2'] # add flow rate

  # add abundance from INFLUX
  id <- findInterval(ALL$time,influx$time)
  if(cruise =="DeepDOM"| cruise == "SCOPE_16" | cruise == "Thompson_9") id <- id +1
  for (phyto in c('pro','syn','pico'))   ALL[,paste0(phyto,'.influx')] <- influx[id,phyto]



  ############### TEST different VC method ###################################
  corr <- 1
  if(cruise == "Thompson_9") corr <- 0.14 # correction factor due to troncated D1 / D2 (based on Seaflow_uncorrected abundance / influx)

  ALL$vc1 <- corr * 1000* 3*ALL[,'fr']*ALL[,'n.opp']/ALL[,'n.evt'] # calculate virtual core v1
  ALL$vc2 <- corr * 1000* 3*ALL[,'fr']*ALL[,'n.opp']/ALL[,'n.evt.'] # calculate virtual core v2
  ALL$vc3 <- corr * 1000* 3*ALL[,'fr']*ALL[,'n.opp.']/ALL[,'n.evt'] # calculate virtual core v3
  ALL$vc4 <- corr * 1000* 3*ALL[,'fr']*ALL[,'n.opp.']/ALL[,'n.evt.'] # calculate virtual core v4



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






############################################
### 3. file-based VC or cruise-based VC  ###
############################################
library(scales)
library(plotrix)

setwd("~/Documents/DATA/Codes/seaflow-virtualcore/2.cruise_calibration/")
DF <- read.csv("SeaflowInflux_comparison.csv")
DF$time <- as.POSIXct(DF$time, tz='GMT')


    s <- 2
   df <- subset(DF, corr==2)
   par(mfrow=c(3,2), pty='s',cex=1.2)
   plot(df[,"pro.influx"], df[,paste0("pro.seaflow.each",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('grey',t), las=1, xlim=c(0,900), ylim=c(0,900)); abline(b=1, a=0, lty=2)
   plot(df[,"pro.influx"], df[,paste0("pro.seaflow.median",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('grey',t), las=1, xlim=c(0,900), ylim=c(0,900)); abline(b=1, a=0, lty=2)
   plot(df[,"syn.influx"], df[,paste0("syn.seaflow.each",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('orange',t), las=1, xlim=c(0,250), ylim=c(0,250)); abline(b=1, a=0, lty=2)
   plot(df[,"syn.influx"], df[,paste0("syn.seaflow.median",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('orange',t),las=1, xlim=c(0,250), ylim=c(0,250)); abline(b=1, a=0, lty=2)
   plot(df[,"pico.influx"], df[,paste0("pico.seaflow.each",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('green',t), las=1, xlim=c(0,100), ylim=c(0,100)); abline(b=1, a=0, lty=2)
   plot(df[,"pico.influx"], df[,paste0("pico.seaflow.median",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('green',t),las=1, xlim=c(0,100), ylim=c(0,100)); abline(b=1, a=0, lty=2)







##################
### 4. BEST VC ###
##################
DF.a <- aggregate(DF, by=list(DF$cruise, DF$corr), mean)
DF.a$cruise <- DF.a$Group.1

# add a color for each cruise
i <- 1
for(cruise in unique(DF$cruise)){
    DF[which(DF$cruise == cruise),'col.cruise'] <- i
    i <- i + 1
  }


VC <- NULL
par(mfrow=c(5,3),cex=1.2, pty='m')
for(c in unique(DF$cruise)){
  print(c)
    df <- subset(DF.a, cruise == c)
    for (phyto in c('pro','syn','pico')){
      #  plot(df$corr, rep(0, nrow(df)), pch=NA ,ylim=c(0,2),xlab='offset',ylab='Diff',main=paste(phyto, c))
      for(s in 1:4){
        # points(df$corr, df[,paste0(phyto, ".diff.each",s)],col=s)
        # points(df$corr, df[,paste0(phyto, ".diff.median",s)],col=s,pch=3)
        each <- mean(df[,paste0(phyto, ".diff.each",s)])
        median <- mean(df[,paste0(phyto, ".diff.median",s)])
        vc <- data.frame(cbind(each, median))
        vc$cruise <- c
        vc$phyto <- phyto
        vc$vc.method <- s
        VC <- rbind(VC, vc)
        }
    }
  }


   met <- aggregate(VC, by=list(VC$vc.method, VC$phyto), FUN=function(x) mean(x, na.rm=T))
   met <- aggregate(VC, by=list(VC$vc.method), FUN=function(x) mean(x, na.rm=T))
   best.vc.method <- unique(c(met[which(met$each == min(met$each)), "vc.method"], met[which(met$median == min(met$median)), "vc.method"]))
   print(paste("best VC method is:",best.vc.method))



######################
### 5. BEST OFFSET ###
######################
s <- best.vc.method
REG <- NULL
for(offset in unique(DF$corr)){
  df <- subset(DF, corr==offset)
  reg.pro <- lm(df[,paste0("pro.seaflow.each",s)] ~ pro.influx, data=df)
  reg.syn <- lm(df[,paste0("syn.seaflow.median",s)] ~ syn.influx, data=df)
  reg.pico <- lm(df[,paste0("pico.seaflow.median",s)] ~ pico.influx, data=df)
  reg <- data.frame(cbind(offset, reg.pro=reg.pro$coefficient[2], reg.syn=reg.syn$coefficient[2], reg.pico=reg.pico$coefficient[2]))
  REG <- rbind(REG, reg)
}
REG$reg.all <- rowMeans(abs(REG[,c("reg.pro", "reg.syn", "reg.pico")]-1))
best.offset <- REG[which(REG$reg.all == min(REG$reg.all)),'offset']
print(paste("best OFFET is:",best.offset))


df <- subset(DF, corr==3)
df2 <- subset(DF, corr==2)
df3 <- subset(DF, corr==4)
par(mfrow=c(3,3), pty='s',cex=1.2)
plot(df[,"pro.influx"], df[,paste0("pro.seaflow.each",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('grey',t), las=1, xlim=c(0,900), ylim=c(0,900)); abline(b=1, a=0, lty=2)
  plot(df2[,"pro.influx"], df2[,paste0("pro.seaflow.each",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('grey',t), las=1, xlim=c(0,900), ylim=c(0,900)); abline(b=1, a=0, lty=2)
  plot(df3[,"pro.influx"], df2[,paste0("pro.seaflow.each",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('grey',t), las=1, xlim=c(0,900), ylim=c(0,900)); abline(b=1, a=0, lty=2)
plot(df[,"syn.influx"], df[,paste0("syn.seaflow.median",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('orange',t),las=1, xlim=c(0,250), ylim=c(0,250)); abline(b=1, a=0, lty=2)
  plot(df2[,"syn.influx"], df2[,paste0("syn.seaflow.median",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('orange',t),las=1, xlim=c(0,250), ylim=c(0,250)); abline(b=1, a=0, lty=2)
  plot(df3[,"syn.influx"], df3[,paste0("syn.seaflow.median",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('orange',t),las=1, xlim=c(0,250), ylim=c(0,250)); abline(b=1, a=0, lty=2)
plot(df[,"pico.influx"], df[,paste0("pico.seaflow.median",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('green',t),las=1, xlim=c(0,100), ylim=c(0,100)); abline(b=1, a=0, lty=2)
  plot(df2[,"pico.influx"], df2[,paste0("pico.seaflow.median",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('green',t),las=1, xlim=c(0,100), ylim=c(0,100)); abline(b=1, a=0, lty=2)
  plot(df3[,"pico.influx"], df3[,paste0("pico.seaflow.median",s)], xlab="Influx", ylab='SeaFlow', pch=21, bg=alpha('green',t),las=1, xlim=c(0,100), ylim=c(0,100)); abline(b=1, a=0, lty=2)




###################
### 6. PLOTTING ###
###################
s <- best.vc.method

png(paste0("SeaFlowInflux-CRUISEcomparison.png"),width=12, height=15, unit='in', res=500)

t <- 0.6
par(mfrow=c(5,3), pty='m',cex=1.2, mar=c(2,3,2,1))
  df1 <- subset(DF, corr==3 | corr==2 | corr==4)

  df2 <- aggregate(df1, by=list(df1$time), FUN=function(x) mean(x, na.rm=T))
    df2b <- aggregate(df1, by=list(df1$time), FUN=function(x) x[1])
    df2$cruise <- df2b$cruise
  df3 <- aggregate(df1, by=list(df1$time), FUN=function(x) sd(x, na.rm=T))
    df3$cruise <- df2b$cruise

  for(c in unique(df2$cruise)){
    df <- subset(df2, cruise == c)
    df.sd <- subset(df3, cruise == c)

    plot(df[,"time"], df[,paste0("pro.influx")], ylab='Pro (cell µL-1)', xlab=NA, pch=1, col=alpha(df$col.cruise,t), ylim=c(min(df[,paste0("pro.influx")], na.rm=T)/2, max(df[,paste0("pro.influx")], na.rm=T)*1.5), las=1)
    plotCI(df[,"time"], df[,paste0("pro.seaflow.each",s)], uiw=df.sd[,paste0("pro.seaflow.each",s)], sfrac=0, pch=21, scol=alpha(df$col.cruise,t), pt.bg=alpha(df$col.cruise,t),add=T)
    plot(df[,"time"], df[,paste0("syn.influx")], ylab='Syn (cell µL-1)', xlab=NA, pch=1, col=alpha(df$col.cruise,t), main=paste(c), ylim=c(min(df[,paste0("syn.influx")], na.rm=T)/2, max(df[,paste0("syn.influx")], na.rm=T)*2), las=1)
    plotCI(df[,"time"], df[,paste0("syn.seaflow.median",s)], uiw=df.sd[,paste0("syn.seaflow.median",s)], sfrac=0, pch=21, scol=alpha(df$col.cruise,t), pt.bg=alpha(df$col.cruise,t),add=T)
    legend('top',legend=c("influx","SeaFlow"), pch=c(1,16),bty='n',cex=0.8)
    if(c == "Thompson_9"){
        plot(df[,"time"], rep(1, nrow(df)), pch=NA)
        legend("center", "NA", bty='n')
      }else{
        plot(df[,"time"], df[,paste0("pico.influx")], ylab='Pico (cell µL-1)', xlab=NA, pch=1, col=alpha(df$col.cruise,t),ylim=c(min(df[,paste0("pico.influx")], na.rm=T)/2, max(df[,paste0("pico.influx")], na.rm=T)*2), las=1)
        plotCI(df[,"time"], df[,paste0("pico.seaflow.median",s)], uiw=df.sd[,paste0("pico.seaflow.median",s)], sfrac=0, pch=21, scol=alpha(df$col.cruise,t), pt.bg=alpha(df$col.cruise,t),add=T)
        }
    print(c)
    print(mean(df[,paste0("pro.seaflow.each",s)], na.rm=T)/mean(df[,paste0("pro.influx")], na.rm=T))
    print(mean(df[,paste0("syn.seaflow.each",s)], na.rm=T)/mean(df[,paste0("syn.influx")], na.rm=T))
    print(mean(df[,paste0("pico.seaflow.median",s)], na.rm=T)/mean(df[,paste0("pico.influx")], na.rm=T))
  }

dev.off()



png(paste0("SeaFlowInflux-ALLcomparison.png"),width=12, height=15, unit='in', res=500)

t <- 0.3
par(mfrow=c(3,2), pty='s',cex=1.2)
  plotCI(df2[,"pro.influx"], df2[,paste0("pro.seaflow.each",s)],df3[,paste0("pro.seaflow.each",s)], sfrac=0, xlab="Influx", ylab='SeaFlow',
        pch=21, scol=df2$col.cruise, pt.bg=alpha(df2$col.cruise,t), main="Prochlorococcus",ylim=c(1,370),xlim=c(1,370), las=1); abline(b=1, a=0, lty=2)
    legend('topleft',legend=unique(df2$cruise),pt.bg=alpha(unique(df2$col.cruise),t), pch=21,bty='n',cex=0.8)
      plotCI(df2[,"pro.influx"], df2[,paste0("pro.seaflow.each",s)],df3[,paste0("pro.seaflow.each",s)], sfrac=0, log='xy', xlab="Influx", ylab='SeaFlow',
        pch=21, scol=df2$col.cruise, pt.bg=alpha(df2$col.cruise,t), main="Prochlorococcus",ylim=c(5,370),xlim=c(5,370), las=1); abline(b=1, a=0, lty=2)
  plotCI(df2[,"syn.influx"], df2[,paste0("syn.seaflow.median",s)],df3[,paste0("syn.seaflow.median",s)], sfrac=0,xlab="Influx", ylab='SeaFlow',
        pch=21, scol=df2$col.cruise, pt.bg=alpha(df2$col.cruise,t), main="Synechococcus",xlim=c(0,150),ylim=c(0,150),las=1); abline(b=1, a=0, lty=2)
      plotCI(df2[,"syn.influx"], df2[,paste0("syn.seaflow.median",s)],df3[,paste0("syn.seaflow.median",s)], sfrac=0, log='xy', xlab="Influx", ylab='SeaFlow',
        pch=21, scol=df2$col.cruise, pt.bg=alpha(df2$col.cruise,t), main="Synechococcus",xlim=c(0.5,150),ylim=c(0.5,150),las=1); abline(b=1, a=0, lty=2)
  plotCI(df2[,"pico.influx"], df2[,paste0("pico.seaflow.median",s)],df3[,paste0("pico.seaflow.median",s)], sfrac=0,xlab="Influx", ylab='SeaFlow',
        pch=21, scol=df2$col.cruise, pt.bg=alpha(df2$col.cruise,t), main="Picoeukaryotes",xlim=c(0,50),ylim=c(0,50),las=1); abline(b=1, a=0, lty=2)
      plotCI(df2[,"pico.influx"], df2[,paste0("pico.seaflow.median",s)],df3[,paste0("pico.seaflow.median",s)], sfrac=0,log='xy', xlab="Influx", ylab='SeaFlow',
        pch=21, scol=df2$col.cruise, pt.bg=alpha(df2$col.cruise,t), main="Picoeukaryotes",xlim=c(0.5,50),ylim=c(0.5,50),las=1); abline(b=1, a=0, lty=2)

dev.off()
