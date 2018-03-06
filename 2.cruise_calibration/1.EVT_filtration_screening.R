library(popcycle)
library(dplyr)
cols <- colorRampPalette(c("blue4", "royalblue4", "deepskyblue3", "seagreen3", "yellow", "orangered2", "darkred"))


######################################
## DOWNLOAD THE RAW DATA using DAT ###
######################################
dat://fa2bb5981465ff583e091071f2f2c5e5f1b118c5219971a2dde107e4dec25e63

#Path to the raw data (DAT)
path.to.data <- "~/Documents/DATA/Codes/seaflow-virtualcore/seaflow-virtualcore-data/"




####################################
### CALCULATE BEADS D1 D2 values ### (only because D1 & D2 values for beads are not recorded in 'stat' table)
####################################

#Path to the Git repository
setwd("~/Documents/DATA/Codes/seaflow-virtualcore/2.cruise_calibration/")




allcruises <- c("SCOPE6", "DeepDOM", "MBARI1","Gradient1")
cruise <- allcruises[4]

for(cruise in allcruises){

 print(cruise)
### Get EVT list
  if(cruise == "SCOPE6" | cruise == "Gradient1") list <- list.files(paste0(path.to.data,cruise,"data"), "00-00$",full.names=T, recursive=T)
  if(cruise == "DeepDOM"| cruise == "MBARI1") list <- list.files(paste0(path.to.data,cruise,"data"), ".evt$",full.names=T, recursive=T)

beads <- NULL
width <- 2500

  for (file in list){

    #file <- list[4]
    print(file)

    evt <- readSeaflow(file,transform=F)

    # Filtering noise
    evt. <- evt[evt$fsc_small > 1 | evt$D1 > 1 | evt$D2 > 1, ]

    # Fltering aligned particles (D1 = D2), with Correction for the difference of sensitivity between D1 and D2
    aligned <- subset(evt., D2 < D1 + width & D1 < D2 + width )



  beads.d1 <- min(subset(aligned, fsc_small > beads.fsc-width/2 & fsc_small < beads.fsc+width/2)[,'D1'])
  beads.d2 <- min(subset(aligned, fsc_small > beads.fsc-width/2 & fsc_small < beads.fsc+width/2)[,'D2'])

  b <- data.frame(cbind(file=basename(file), beads.d1, beads.d2))
  beads <- rbind(beads, b)

  }

  write.csv(beads, paste0(cruise,"data/beadsD1D2.csv"), quote=F, row.names=F)
}






#################################
### SET FILTRATION PARAMETERS ###
#################################
gate <- FALSE

allcruises <- c("SCOPE6", "DeepDOM", "MBARI1","Gradient1")

cruise <- allcruises[2]

for(cruise in allcruises){

print(cruise)
### Get EVT list
if(cruise == "SCOPE6" | cruise == "Gradient1") list <- list.files(paste0(path.to.data,cruise,"data"), "00-00$",full.names=T, recursive=T)
if(cruise == "DeepDOM"| cruise == "MBARI1") list <- list.files(paste0(path.to.data,cruise,"data"), ".evt$",full.names=T, recursive=T)



  if(cruise == "SCOPE6") inst <- 740
  if(cruise == "DeepDOM" | cruise == "MBARI1") inst <- 989
  if(cruise == "Gradient1") inst <- 751

  # SLOPES
  slopes <- read.csv("~/Documents/DATA/Codes/seaflow-virtualcore/1.bead_calibration/seaflow_filter_MERGED_slopes.csv")
    notch.small.D1 <- slopes[slopes$ins== inst,'notch.small.D1']
    notch.small.D2 <- slopes[slopes$ins== inst,'notch.small.D2']
    notch.large.D1 <- slopes[slopes$ins== inst,'notch.large.D1']
    notch.large.D2 <- slopes[slopes$ins== inst,'notch.large.D2']

  # BEADS
    if(cruise == "SCOPE6"){
      stat <- read.csv(paste0(cruise,"data/stats.csv"))
      stat$time <- as.POSIXct(stat$time, format = "%FT%T", tz = "GMT")
      beads <- subset(stat, pop == 'beads')
      beads$fsc_small <- (log10(beads$fsc_small)/3.5)*2^16# LOG TRANSFORM DATA

      beads.fsc <- 58000
      plot(beads$time, beads$fsc_small,ylim=c(0,2^16))
      abline(h=beads.fsc, col=1)
      print(beads.fsc)
    }

    if(cruise == "MBARI1"){
      stat <- read.csv(paste0(cruise,"data/stats.csv"))
      stat$time <- as.POSIXct(stat$time,  tz = "GMT")
      beads <- subset(stat, pop == 'beads')
      beads.fsc <- median(beads$fsc_small)
      plot(beads$time, beads$fsc_small,ylim=c(0,2^16))
      abline(h=beads.fsc, col=2)
      print(beads.fsc)
    }

    if(cruise == "DeepDOM" | cruise == "Gradient1"){
      stat <- read.csv(paste0(cruise,"data/stats.csv"))
      stat$time <- as.POSIXct(stat$time, format = "%FT%T", tz = "GMT")
      beads <- subset(stat, pop == 'beads')
      beads$fsc_small <- (log10(beads$fsc_small)/3.5)*2^16# LOG TRANSFORM DATA

      beads.fsc <- median(beads$fsc_small)
      plot(beads$time, beads$fsc_small,ylim=c(0,2^16))
      abline(h=beads.fsc, col=2)
      print(beads.fsc)
    }


    # CALCULATE BEADS D1 D2 values (only because D1 & D2 values for beads are not recorded in 'stat' table)
    beads <- read.csv(paste0(cruise,"data/beadsD1D2.csv"))
    beads.d1 <- median(beads$beads.d1)
    beads.d2 <- median(beads$beads.d2)

    # plot(beads$beads.d1, ylim=c(0,2^16)); abline(h=beads.d1); abline(h=beads.d2, col=2)
    #   points(beads$beads.d2,col=2)

      correction.D1 <- round(beads.d1 - notch.small.D1 * beads.fsc)
      correction.D2 <- round(beads.d2 - notch.small.D2 * beads.fsc)


width <- 2500
ALL <- NULL

for (file in list){

  #file <- list[18]
  print(file)

  evt <- readSeaflow(file,transform=F)

  # Filtering noise
  evt. <- evt[evt$fsc_small > 1 | evt$D1 > 1 | evt$D2 > 1, ]

  # Fltering aligned particles (D1 = D2), with Correction for the difference of sensitivity between D1 and D2
  aligned <- subset(evt., D2 < D1 + width & D1 < D2 + width)


    if(nrow(aligned)> 100000){aligned. <- sample_n(aligned, 100000)
      } else aligned. <- aligned

    # if(nrow(evt.)> 100000){df <- sample_n(evt., 100000)
    #   } else df <- evt.
    # par(mfrow=c(2,2))
    # plot.cytogram(df, "D1", "D2"); abline(b=notch, a=origin)
    # plot.cytogram(aligned., "D1", "D2"); abline(b=1, a=origin)
    # plot.cytogram(df, "fsc_small", "D1"); abline(b=notch.small.D1, a=offset.small.D1,col=2); abline(b=notch.large.D1, a=offset.large.D1,col=3)
    # plot.cytogram(df, "fsc_small", "D2"); abline(b=notch.small.D2, a=offset.small.D2,col=2); abline(b=notch.large.D2, a=offset.large.D2,col=3)


### 1. intial screening
# screening <- c(correction, -20000,17500,-15000,-12500,-10000,-7500,-5000,-2500,0,2500,5000,7500,10000, 12500, 15000)

### 2. offet for small particles set to 0 (and corretced offset of large particles )
# screening <- c(-10000,-7500,-5000,-2500,0,2500,5000,7500,10000)

### 3. fined tune screening (offet for small particles set to 0, same as 2.):
screening <- seq(-2500,2500, by=500)

  for(corr in screening){

    print(paste("corr=", corr))
    #corr <- 0
  ### OFFSET intial screening
  # offset.small.D1 <- round(beads.d1 - notch.small.D1 * beads.fsc) + corr
  # offset.small.D2 <- round(beads.d2 - notch.small.D2 * beads.fsc) + corr
  # offset.large.D1 <- round(beads.d1 - notch.large.D1 * beads.fsc) + corr
  # offset.large.D2 <- round(beads.d2 - notch.large.D2 * beads.fsc) + corr
  ### Corrected OFFSET
  offset.small.D1 <- round(beads.d1 - notch.small.D1 * beads.fsc - correction.D1) + corr
  offset.small.D2 <- round(beads.d2 - notch.small.D2 * beads.fsc - correction.D2) + corr
  offset.large.D1 <- round(beads.d1 - notch.large.D1 * beads.fsc - correction.D1) + corr
  offset.large.D2 <- round(beads.d2 - notch.large.D2 * beads.fsc - correction.D2) + corr

opp <- subset(aligned, D1 <= fsc_small*notch.small.D1 + offset.small.D1 & D2 <= fsc_small*notch.small.D2 + offset.small.D2 |
    D1  <= fsc_small*notch.large.D1 + offset.large.D1 & D2 <= fsc_small*notch.large.D2 + offset.large.D2)
    if(nrow(opp) < 10) next

opp$pop <- 0

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




# png(paste0(cruise,"data/", basename(file),"-corr",corr,".png"),width=9, height=12, unit='in', res=100)

# if(nrow(opp)> 100000){opp. <- sample_n(opp, 100000)
#    } else opp. <- opp
# par(mfrow=c(2,2))
# plot.cytogram(aligned. ,  "fsc_small", "D1"); abline(b=notch.small.D1, a=offset.small.D1,col=2); abline(b=notch.large.D1, a=offset.large.D1,col=3)
# plot.cytogram(aligned. , "fsc_small", "D2"); abline(b=notch.small.D2, a=offset.small.D2,col=2); abline(b=notch.large.D2, a=offset.large.D2,col=3)
# plot.vct.cytogram(opp. , "fsc_small", "pe")
# plot.vct.cytogram(opp. , "fsc_small", "chl_small")

# dev.off()

n.opp <- nrow(opp)
n.evt <- nrow(evt)
n.evt. <- nrow(evt.)
beads <- nrow(b)
syn <- nrow(s)
pro <- nrow(p)
pico <- nrow(l)
filename <- basename(file)
if(cruise == 'MBARI1' | cruise == 'DeepDOM') filename <- paste(basename(dirname(file)),basename(file),sep="/")

all <- data.frame(cbind(file=filename, n.opp, n.evt, n.evt., beads, pro, syn, pico, origin, offset.small.D1, offset.large.D1,offset.small.D2, offset.large.D2, width, corr))

ALL <- rbind(ALL, all)

    }
  }

  write.csv(ALL, paste0(cruise,"data/seaflow-summary3.csv"), quote=F, row.names=F)


}















###############
### SUMMARY ###
###############
library(lattice)
library(lmodel2)

s <- 3#1,2,3 ...

setwd("~/Documents/DATA/Codes/seaflow-virtualcore/2.cruise_calibration/")

allcruises <- c("SCOPE6", "DeepDOM", "MBARI1","Gradient1")
cruise <- allcruises[4]

summary <- NULL
for(cruise in allcruises){

print(cruise)

### SFL + INFLUX
if(cruise == 'DeepDOM'){sfl <- read.csv(paste0(cruise,"data/sfl.csv"))
                        sfl$date <- as.POSIXct(sfl$date, format = "%FT%T", tz = "GMT")
                      influx <- read.delim(paste0(cruise,"data/surface_samples_metadata.tab"))
                        influx$time <- as.POSIXct(influx$Date.Time, format="%m/%d/%y %H:%M", tz="GMT")
                        }
if(cruise == 'Gradient1'){sfl <- read.csv(paste0(cruise,"data/sfl.csv"))
                            sfl$date <- as.POSIXct(sfl$date, format = "%FT%T", tz = "GMT")
                            sfl$file <- basename(as.character(sfl$file))
                          influx <- read.csv(paste0(cruise,"data/all_station_curated.csv"))
                            influx$time_HST <- as.POSIXct(influx$time_HST, format="%m/%d/%y %H:%M", tz="HST")
                            influx <- influx[order(influx$time_HST), ]
                            influx <- subset(influx, Depth_m < 16)
                            pro <- subset(influx, population == 'prochloro')
                            syn <- subset(influx, population == 'synecho')
                            picoeuk <- subset(influx, population == 'picoeuk')
                          influx <- data.frame(cbind(time=pro$time_HST, pro=pro$abundance, syn=syn$abundance, pico=picoeuk$abundance))
                            influx$time <- as.POSIXct(influx$time, origin='1970-01-01',tz="GMT")
                              }
if(cruise == 'MBARI1'){sfl <- read.delim(paste0(cruise,"data/cruise.sfl"))
                        sfl$date <- as.POSIXct(sfl$DATE, format = "%FT%T", tz = "GMT")
                        sfl$file <- sfl[,"FILE"]
                        sfl$flow_rate <- sfl[,"FLOW.RATE"]
                      influx <- read.csv(paste0(cruise,"data/mbari_summary.csv"))
                        id <- match(influx[,'file'],sfl[,"FILE"])
                        influx$time <- sfl[id, 'date']
                        influx$pro <- influx$pro / 1000
                        influx$syn <- influx$syn / 1000
                        influx$pico <- influx$picoeuks / 1000
                          }
if(cruise == 'SCOPE6'){sfl <- read.csv(paste0(cruise,"data/sfl.csv"))
                        sfl$date <- as.POSIXct(sfl$date, format = "%FT%T", tz = "GMT")
                        sfl$flow_rate <- median(sfl$flow_rate, na.rm=T)
                        sfl$file <- basename(as.character(sfl$file))
                      influx <- read.csv(paste0(cruise,"data/influx.csv"))
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
ALL <- read.csv(paste0(cruise,"data/seaflow-summary",s,".csv"))
if(cruise == "Gradient1") ALL <- ALL[!(ALL$file =="2016-04-26T15-07-38-00-00"),]
if(cruise == "DeepDOM") ALL <- ALL[!(ALL$file =="2013_094/321.evt" | ALL$file=="2013_124/409.evt"),]

  id2 <- match(ALL[,'file'],as.character(sfl[,"file"]))
  ALL$time <- sfl[id2, 'date']
  ALL$fr <- sfl[id2, 'flow_rate']
  ALL <- ALL[order(ALL$time), ]
  ALL$vc <- 3*ALL[,'fr']*ALL[,'n.opp']/ALL[,'n.evt.']
  ALL[which(ALL$pro == 120),'pro'] <- 1


# plot(ALL$corr, ALL$vc,pch=2, log='y')

id <- findInterval(ALL$time,influx$time)
if(cruise =="DeepDOM"| cruise == "Gradient1") id <- id +1
ALL[,'pro.influx'] <- influx[id,'pro']
ALL[,'syn.influx'] <- influx[id,'syn']
ALL[,'pico.influx'] <- influx[id,'pico']

  # if(cruise == "SCOPE6")  VC <- 0.136
  # if(cruise == "DeepDOM" | cruise == "MBARI1") VC <- 0.149
  # if(cruise == "Gradient1") VC <- 0.143


  # ### PLOTING ALL
  # par(mfrow=c(3,1),cex=1.2, pty='m')
  # plot(ALL$time, ALL[,'pro']/ALL[,'vc'], col=1, main="pro", ylab="Abundance", xlab=NA)
  #    points(ALL$time,ALL[,'pro.influx'],pch=1, cex=2, col='darkgrey',type='o')
  # plot(ALL$time, ALL[,'syn']/ALL[,'vc'], col=1, main="syn", ylab="Abundance", xlab=NA)
  #   points(ALL$time,ALL[,'syn.influx'],pch=1, cex=2, col='darkgrey',type='o')
  # plot(ALL$time, ALL[,'pico']/ALL[,'vc'], col=1, main="pico", ylab="Abundance", xlab=NA)
  #   points(ALL$time,ALL[,'pico.influx'],pch=1, cex=2, col='darkgrey',type='o')


### LOOK FOR THE BEST OFFSET
results <- results2 <- NULL
for(offset in unique(ALL$corr)){
  df <- subset(ALL,corr == offset, sep=',')
    df[,'pro.seaflow'] <- df[,'pro']/df[,'vc']
    df[,'pro.seaflow2'] <- df[,'pro']/median(df[,'vc'])
    df[,'syn.seaflow'] <- df[,'syn']/df[,'vc']
    df[,'syn.seaflow2'] <- df[,'syn']/median(df[,'vc'])
    df[,'pico.seaflow'] <- df[,'pico']/df[,'vc']
    df[,'pico.seaflow2'] <- df[,'pico']/median(df[,'vc'])

    res <- data.frame(offset, vc=c('each','median'), diff.pro=rbind(sum(abs(df$pro.influx-df$pro.seaflow), na.rm=T), sum(abs(df$pro.influx-df$pro.seaflow2), na.rm=T)),
                                                    diff.syn=rbind(sum(abs(df$syn.influx-df$syn.seaflow), na.rm=T), sum(abs(df$syn.influx-df$syn.seaflow2), na.rm=T)),
                                                    diff.pico=rbind(sum(abs(df$pico.influx-df$pico.seaflow), na.rm=T), sum(abs(df$pico.influx-df$pico.seaflow2), na.rm=T)))
      results <- rbind(results, res)
      }
print(results)

par(mfrow=c(3,1),cex=1.2, mar=c(4,4,2,1), pty='m')
plot(results$offset, results$diff.pro, col=results$vc)
plot(results$offset, results$diff.syn, col=results$vc)
plot(results$offset, results$diff.pico, col=results$vc)





    if(s ==1){
      if(cruise == "SCOPE6"){ id.pro <- 1 ; id.syn <- id.pico <- 2 }
      if(cruise == "DeepDOM"){ id.pro <- 1 ; id.syn <- id.pico <- 2 }
      if(cruise == "MBARI1"){ id.pro <- 1 ; id.syn <- id.pico <- 2 }
      if(cruise == "Gradient1"){ id.pro <-  1 ; id.syn <- id.pico <- 2 }
    }

    if(s ==2){
      if(cruise == "SCOPE6"){ id.pro  <- which(results$offset == 0)[1]; id.syn <- id.pico <- id.pro+1}
      if(cruise == "DeepDOM"){ id.pro <- which(results$offset == 0)[1]; id.syn <- id.pico <- id.pro+1}
      if(cruise == "MBARI1"){ id.pro  <- which(results$offset == 0)[1]; id.syn <- id.pico <- id.pro+1 }
      if(cruise == "Gradient1"){ id.pro <- which(results$offset == 0)[1]; id.syn <-id.pico <- id.pro+1 }
  }

    if(s==3){
      if(cruise == "SCOPE6"){ id.pro  <- which(results$offset == 0)[1] ; id.syn <- id.pico <- id.pro+1}
      if(cruise == "DeepDOM"){ id.pro <- which(results$offset == 0)[1]; id.syn <- id.pico <- id.pro+1}
      if(cruise == "MBARI1"){ id.pro  <- which(results$offset == 0)[1]; id.syn <- id.pico <- id.pro+1 }
      if(cruise == "Gradient1"){ id.pro <- which(results$offset == 0)[1]; id.syn <- id.pico <- id.pro+1 }
    }


    df.pro <- subset(ALL, corr == results[id.pro,'offset'])
    df.syn <- subset(ALL, corr == results[id.syn,'offset'])
    df.pico <- subset(ALL, corr == results[id.pico,'offset'])

    DF <- data.frame(cbind(cruise= cruise, pro.influx=df.pro[,"pro.influx"], syn.influx=df.syn[,"syn.influx"], pico.influx=df.pico[,"pico.influx"],
                                           pro.seaflow=df.pro[,"pro"]/df.pro[,'vc'], syn.seaflow=df.syn[,"syn"]/median(df.syn[,'vc']), pico.seaflow=df.pico[,"pico"]/median(df.pico[,'vc'])))
    summary <- rbind(summary, DF)



png(paste0(cruise, "-SeaFlowInflux-comp",s,".png"),width=12, height=12, unit='in', res=100)

    par(mfrow=c(3,2),cex=1.2, mar=c(4,4,2,1))


    df <- subset(ALL, corr == results[id.pro,'offset'])
    par(pty='m')
    lim <- 300
    plot(df$time, df[,'pro']/df[,'vc'], main=paste("pro",results[id.pro,'offset'], results[id.pro,'vc']), ylim=c(0,lim), ylab="Abundance", xlab=NA,type='o')
      points(df$time, df[,'pro']/median(df[,'vc']), type='o',col=2)
      points(df$time, df[,'pro.influx'],pch=16, cex=2, col='darkgrey',type='o')

    par(pty='s')
    plot(df[,'pro']/df[,'vc'],df[,'pro.influx'],ylab="Influx", xlab="SeaFlow",ylim=c(10,lim),xlim=c(0,lim)); abline(b=1,a=0)

    df <- subset(ALL, corr == results[id.syn,'offset'])
    par(pty='m')
    lim <- 150
    plot(df$time, df[,'syn']/df[,'vc'], main=paste("syn",results[id.syn,'offset'], results[id.syn,'vc']), ylab="Abundance", xlab=NA,type='o')
      points(df$time, df[,'syn']/median(df[,'vc']), type='o',col=2)
      points(df$time, df[,'syn.influx'],pch=16, cex=2, col='darkgrey',type='o')

    par(pty='s')
    plot(df[,'syn']/median(df[,'vc']),df[,'syn.influx'],ylab="Influx", xlab="SeaFlow"); abline(b=1,a=0)

    df <- subset(ALL, corr == results[id.pico,'offset'])
    par(pty='m')
    lim <- 40
    plot(df$time, df[,'pico']/df[,'vc'], main=paste("pico",results[id.pico,'offset'], results[id.pico,'vc']), ylab="Abundance", xlab=NA,type='o')
      points(df$time, df[,'pico']/median(df[,'vc']), type='o',col=2)
      points(df$time, df[,'pico.influx'],pch=16, cex=2, col='darkgrey',type='o')

    par(pty='s')
    plot(df[,'pico']/median(df[,'vc']),df[,'pico.influx'], ylab="Influx", xlab="SeaFlow"); abline(b=1,a=0)

  dev.off()



  # results <- results[order(results$offset),]
  # results <- results[order(results$vc),]
  # par(mfrow=c(1,1))
  # plot(results$offset, results$diff.pro, log='y', ylim=c(min(results$diff.pico),max(results$diff.pro)),pch=as.numeric(results$vc),type='o')
  # points(results$offset, results$diff.syn, col=2, pch=as.numeric(results$vc),type='o')
  # points(results$offset, results$diff.pico, col=3,pch=as.numeric(results$vc),type='o')




}

write.csv(summary,paste0("SeaFlow-Influx_compALL",s,".csv"), quote=F, row.names=F)






############
### BEST ###
############
library(lmodel2)
setwd("~/Documents/DATA/Codes/seaflow-virtualcore/2.cruise_calibration/")
summary <- read.csv(paste0( "SeaFlow-Influx_compALL",s,".csv"))

png(paste0("SeaFlowInflux-compALL",s,".png"),width=6, height=12, unit='in', res=100)
  par(mfrow=c(3,1),cex=1.2, mar=c(4,4,2,1),pty='s')
  plot(summary[,'pro.seaflow'],summary[,'pro.influx'],asp=1,ylab="Influx", xlab="SeaFlow",col=summary[,'cruise'], pch=16, ylim=c(50,500), xlim=c(50,500) ,log='xy'); abline(b=1,a=0,lwd=2,col='grey')
  plot(summary[,'syn.seaflow'],summary[,'syn.influx'],asp=1,ylab="Influx", xlab="SeaFlow",col=summary[,'cruise'], pch=16, ylim=c(0.5,200), xlim=c(0.5,200), log='xy'); abline(b=1,a=0,lwd=2,col='grey')
  plot(summary[,'pico.seaflow'],summary[,'pico.influx'],asp=1,ylab="Influx", xlab="SeaFlow",col=summary[,'cruise'], pch=16, ylim=c(0.5,100), xlim=c(0.5,100), log='xy'); abline(b=1,a=0,lwd=2,col='grey')
dev.off()
