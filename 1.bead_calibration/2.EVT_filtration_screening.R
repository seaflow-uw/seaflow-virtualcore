library(popcycle)
library(dplyr)


######################################
## DOWNLOAD THE RAW DATA using DAT ###
######################################
#dat://fa2bb5981465ff583e091071f2f2c5e5f1b118c5219971a2dde107e4dec25e63

#Path to the raw data (DAT)
path.to.data <- "~/Documents/DATA/Codes/seaflow-virtualcore/seaflow-virtualcore-data/"



########################
### SEAFLOW ANALYSIS ###
########################

#Path to the Git repository
setwd("~/Documents/DATA/Codes/seaflow-virtualcore/")

for(inst in c(989,751,740)){
  print(paste("instrument=",inst))

  list <- list.files(path=paste0(path.to.data,inst,"_caldata"),pattern=".evt$", recursive=T, full.names=T)
  if(inst == 740){
    list <- list[-14]
    pe.lim.o <- 50000
    fsc.lim.o <- 45000
    pe.lim.c <- 35000
    fsc.lim.c <- 40000
    }

  if(inst == 751){
    list <- list[-14]
    pe.lim.o <- 55000
    fsc.lim.o <- 50000
    pe.lim.c <- 10000
    fsc.lim.c <- 37500
  }

  if(inst == 989){
    list <- list[-c(15:16)]
    pe.lim.o <- 50000
    fsc.lim.o <- 45000
    pe.lim.c <- 35000
    fsc.lim.c <- 35000
  }


  slopes <- read.csv("1.bead_calibration/seaflow_filter_MERGED_slopes.csv")
    notch.small.D1 <- slopes[slopes$ins== inst,'notch.small.D1']
    notch.small.D2 <- slopes[slopes$ins== inst,'notch.small.D2']
    notch.large.D1 <- slopes[slopes$ins== inst,'notch.large.D1']
    notch.large.D2 <- slopes[slopes$ins== inst,'notch.large.D2']


  ## BEADS
    ALL <- read.csv(paste0("1.bead_calibration/",inst,"-summary.csv"))
      ALL$pmt <- 'coastal'
      id <- which((grepl("oligo", ALL$file))==T)
      ALL[id, "pmt"] <- 'oligo'
      ALL$size <- as.numeric(sub("um.evt","",t(data.frame(list(strsplit(as.character(ALL$file),"_"))))[,5], ))
      ALL <- ALL[order(ALL$size),]

     ref.o <- subset(ALL, pmt == 'oligo'  & size == 1)
     ref.c <- subset(ALL, pmt == 'coastal'  & size == 1)
  if(grepl("oligo", list[i])){ beads.fsc <- ref.o$fsc.med
                                beads.d1 <- ref.o$d1.med
                                beads.d2 <- ref.o$d2.med
                    }else{      beads.fsc <- ref.c$fsc.med
                                beads.d1 <- ref.c$d1.med
                                beads.d2 <- ref.c$d2.med
                    }

  width <- 3000

  screening <- c(0,2500,5000,7500,10000, 1, 2, 3, 4)


  ALL <- NULL
  for(corr in screening){

  print(paste("corr=", corr))

  for (i in 1:length(list)){
  #i <- 10

  print(list[i])

  evt <- readSeaflow(list[i],transform=F)

    # Filtering particles detected by fsc_small
    evt. <- evt[evt$fsc_small > 1 | evt$D1 > 1 | evt$D2 > 1, ]
    evt.2 <- evt.[evt.$fsc_small > 6000, ]

    # Fltering aligned particles (D1 = D2)
    aligned <- subset(evt.2, D2 < D1 + width & D1 < D2 + width )


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

  if(nrow(opp) > 10){

    b <- opp[which(opp$fsc_small > 6000),] #beads
    if(basename(list[i]) == "2017-02-21_740_oligo_5.7um.evt" | basename(list[i]) == "2017-02-28_989_oligo_5.7um.evt") b <- opp[which(opp$chl_small > 50000),]
    beads <- nrow(b)
    opp.evt.ratio <- round(nrow(opp)/nrow(evt.),3)
    opp.evt.ratio2 <- round(nrow(opp)/nrow(evt.2),3)

      if(nrow(opp)> 100000){ opp. <- sample_n(opp, 100000)
      } else opp. <- opp

          if(corr == 0 | corr == 2 | corr == 3 | corr == 4){
          png(paste0("1.bead_calibration/",inst, "_caldata/", basename(list[i]),"corr",corr,".png"),width=9, height=12, unit='in', res=100)
            par(mfrow=c(2,2))
            plot.cytogram(opp. ,  "fsc_small", "D1"); abline(b=notch.small.D1, a=offset.small.D1,col=2); abline(b=notch.large.D1, a=offset.large.D1,col=3)#; points(beads$fsc_small, beads$D1, col=2, pch=16)
            plot.cytogram(opp. , "fsc_small", "D2"); abline(b=notch.small.D2, a=offset.small.D2,col=2); abline(b=notch.large.D2, a=offset.large.D2,col=3)#; points(beads$fsc_small, beads$D2, col=2, pch=16)
            plot.cytogram(opp. , "fsc_small", "pe"); # points(beads$fsc_small, beads$pe, col=2, pch=1)
            plot.cytogram(opp. , "fsc_small", "chl_small"); #points(beads$fsc_small, beads$chl_small, col=2, pch=1)
          dev.off()
          }

    }else{ beads <- 0
          opp.evt.ratio <- opp.evt.ratio2 <- NA
          }

  all <- data.frame(cbind(file=as.character(list[i]), beads, opp.evt.ratio, opp.evt.ratio2, width=width, offset=corr))

  ALL <- rbind(ALL, all)

  }
  }

  write.csv(ALL, paste0("1.bead_calibration/",inst,"-beads_count.csv"), quote=F, row.names=F)

  }











  ###########
  ### REF ###
  ###########
  for(ins in c(989,751,740)){

  ref <- read.csv("1.bead_calibration/beadcal_data_FR.csv")
  ref. <- subset(ref, inst == ins)
  ref.$col <- 1
  id <- which((grepl("oligo", ref.$file))==T)
  ref.[id, 'col'] <- 2


  ddr <- read.csv("detectable_region.csv")
  VC <- mean(drr[which(ddr$seaflow_serial == ins), "detectable_region_ratio"])

  fr <- read.csv("1.bead_calibration/beadcalruntimes_flowrates.csv")
  fr$filename <- paste0(ins,"_caldata/",fr$FILE.NAME,".evt")
  fr$beads.size <- as.numeric(sub("um.evt","",t(data.frame(list(strsplit(fr$file,"_"))))[,5], ))

  png(paste0("1.bead_calibration/",ins,"-Influx-vs-SeaFlow.png"),width=12, height=9, unit='in', res=100)
  par(mfrow=c(1,2), pty='s')

  #for(corr in c(2500,5000,7500,10000)){
  corr <- 7500
    ALL <- read.csv(paste0("1.bead_calibration/",ins,"-beads_count.csv"))
    #ALL$beads.size <- as.numeric(sub("um.evt","",t(data.frame(list(strsplit(as.character(ALL$file),"_"))))[,5], ))
        # plot(ALL$beads.size, ALL$opp.evt.ratio, col=round(ALL$offset/2000))
    ALL <- subset(ALL, offset == corr)
    merge <- merge(fr, ALL, by.x="filename", by.y="file")
    merge$abundance <- round(merge$beads/(merge$volume..ml.*VC*merge$opp.evt.ratio))/1000
    merge$abundance2 <- round(merge$beads/(merge$volume..ml.*VC*merge$opp.evt.ratio2))/1000
    merge$col <- 1
    id <- which((grepl("oligo", merge$filename))==T)
    merge[id, 'col'] <- 2

  df <- merge(merge, ref., by.y=c("Bead.Size..um.","col"), by.x=c("beads.size",'col'))
  df <- subset(df, beads.size > 0.3)
  plot(df$abundance.y, df$abundance2, pch=NA,asp=1, main=paste(ins, "offset =",corr), log='xy', ylim=c(1,500), xlim=c(1,500))
  text(df$abundance.y, df$abundance2, labels=df$beads.size, col=df$col+2)
  abline(0,1)
  # plot(df$abundance.y, df$abundance.x, pch=NA,asp=1, main=paste(ins, "offset =",corr), log='xy', ylim=c(1,500), xlim=c(1,500))
  # text(df$abundance.y, df$abundance.x, labels=df$beads.size, col=df$col)
  # abline(0,1)

      plot(ref.$Bead.Size..um., ref.$abundance, col=ref.$col, cex=2, lwd=2, pch=3, main=paste(ins, "offset =",corr), ylab="Abundance", log='y',xlab="Size", ylim=c(1,500))
      # points(merge$beads.size, merge$abundance, cex=2, col=merge$col)
      points(merge$beads.size, merge$abundance2, cex=2, col=merge$col+2)
      legend("top", c("coastal","oligo"), text.col=c(1,2), bty='n')

dev.off()

}



###############
#### SLOPES ###
###############
setwd("~/Documents/DATA/Codes/seaflow-virtualcore")

png("1.bead_calibration/ALL_MERGED_Slopes.png",width=12, height=15, unit='in', res=100)
par(mfrow=c(3,2),cex=1.2)

SLOPES <- NULL

for(ins in c(989,751,740)){
print(ins)

      #ins<- "751"
    ALL <- read.csv(paste0("1.bead_calibration/",ins,"-summary.csv"))
      ALL$pmt <- 'coastal'
      id <- which((grepl("oligo", ALL$file))==T)
      ALL[id, "pmt"] <- 'oligo'
      ALL$beads.size <- as.numeric(sub("um.evt","",t(data.frame(list(strsplit(as.character(ALL$file),"_"))))[,5], ))
      ALL <- ALL[order(ALL$beads.size),]

       df.o <- subset(ALL, pmt == 'oligo' & beads > 0 & beads.size > 0.18)
       ref.o <- subset(df.o, beads.size ==1.00)
       df.c <- subset(ALL, pmt == 'coastal' & beads > 0 & beads.size > 0.18)
       ref.c <- subset(df.c, beads.size ==1.00)

       if(ins == 740){
         comp.fsc <- mean(df.o$fsc.med[-1] - df.c$fsc.med)# [-1] to remove data from 0.3 µm beads not analyzed with instrument 740 with Coastal settings
         comp.d1 <- mean(df.o$d1.med[-1] - df.c$d1.med)
         comp.d2 <- mean(df.o$d2.med[-1] - df.c$d2.med)
         }else{
         comp.fsc <- mean(df.o$fsc.med - df.c$fsc.med)
         comp.d1 <- mean(df.o$d1.med - df.c$d1.med)
         comp.d2 <- mean(df.o$d2.med - df.c$d2.med)
       }

    DF <- data.frame(cbind(fsc=c(df.o$fsc.med,df.c$fsc.med+comp.fsc), D1=c(df.o$d1.med,df.c$d1.med+comp.d1),D2=c(df.o$d2.med,df.c$d2.med+comp.d2), size=c(df.o$beads.size,df.c$beads.size)))
    DF <- DF[order(DF$size),]


    slope <- NULL
    for(size in c("small","large")){
      print(size)
      if(size == "small" ) sm <- subset(DF, size <= 1)
      if(size == "large" ) sm <- subset(DF, size >= 1)

      for(param in c("D1","D2")){
        reg <- data.frame(cbind(fsc=sm$fsc, predict(lm(sm[,param] ~ fsc, data=sm), newdata=data.frame(fsc=sm$fsc),interval='predict')))

          for(q in c("fit","lwr","upr")){
            pred <- t(rev(round(lm(reg[,q] ~ fsc, data=reg)$coefficient,3)))
              if(q == "fit") colnames(pred) <- c(paste0("notch.",size,".",param),paste0("intersect.",size,".",param))
              if(q == "lwr") colnames(pred) <- c(paste0("notch.",size,".",param,"_2.5"),paste0("intersect.",size,".",param,"_2.5"))
              if(q == "upr") colnames(pred) <- c(paste0("notch.",size,".",param,"_97.5"),paste0("intersect.",size,".",param,"_97.5"))
            slope <- data.frame(cbind(slope, pred))
          }
        }
      }



      for(param in c("D1","D2")){
        col <- c("red3",rep('grey',2),"red3",rep('grey',2))
        plot(0,0,pch=NA, xlim=c(0,62000), ylim=c(0,62000),main=paste(ins), ylab=paste(param), xlab="FSC")

          s <- grep(param, colnames(slope))
            c <- 1
            for(i in s[seq(1,length(s),by=2)]){
                abline(b=slope[,i],a=slope[,i+1], lwd=2, col=col[c])
                c <- c+1}

                points(DF$fsc, DF[,param], pch=21, bg=rep(alpha(viridis(7),0.5),each=2), cex=2)
                legend("topleft",legend=paste(unique(DF$size), "µm beads"), cex=1, pch=21, bty='n',pt.bg=alpha(viridis(7),0.5))

              }


  slope <- data.frame(cbind(ins=ins, slope))
  SLOPES <- rbind(SLOPES, slope)

}

dev.off()

write.csv(SLOPES, "1.bead_calibration/seaflow_filter_slopes.csv", quote=F, row.names=F)
write.csv(DF, paste0("1.bead_calibration/",inst,"-beadscoordinates.csv"), quote=F, row.names=F)
