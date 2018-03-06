library(flowCore)
library(splancs)
library(plotrix)
library(caroline)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))


plot.cytogram <- function (evtopp, para.x = "FSC.small.stuff", para.y = "X692.40.small.stuff", ...)
{
    cols <- colorRampPalette(c("blue4", "royalblue4", "deepskyblue3",
        "seagreen3", "yellow", "orangered2", "darkred"))
    par(pty = "s")
    plot(evtopp[, c(para.x, para.y)], pch = 16, cex = 0.3,col = densCols(log10(evtopp[, c(para.x, para.y)]),
          colramp = cols), xlim = c(1, 10^4), ylim = c(1, 10^4), log='xy',...)
          }

plot.vct.cytogram <- function (opp, para.x = "fsc_small", para.y = "chl_small", ...)
          {
        plot(opp[, c(para.x, para.y)], pch = 16, cex = 0.3, col = as.numeric(as.factor(opp$pop)),
                    xlim = c(1, 10^4), ylim = c(1, 10^4), log='xy', ...)
              legend("topleft", legend = (unique(opp$pop)), col = unique(as.numeric(as.factor(opp$pop))),
                  pch = 16, pt.cex = 0.6, bty = "n")
                    }


#################
## BATCH FILES ##
#################

setwd("~/Documents/DATA/Codes/seaflow-virtualcore")

savepath <- "1.bead_calibration/Influx_bead_data"


file.list <- dir("seaflow-virtualcore-data/Influx_bead_data", pattern = ".fcs$", recursive=F, full.names=T)

summary.table <- NULL
draw.gate <- TRUE

for (file in file.list) {
    print(paste("processing file:",file))


#file <- file.list[1]
###############
## read FCS ###
###############

fcs <- read.FCS(file, transformation=T, emptyValue=F)
opp <- tab2df(exprs(fcs))
#opp$pop <- 0


# par(mfrow=c(2,2))
# plot.cytogram(opp, "FSC.small.stuff","X692.40.small.stuff")
# plot.cytogram(opp, "FSC.small.stuff","X580.30")
# plot.cytogram(opp, "X692.40.small.stuff","X580.30")
# plot.cytogram(opp, "SSC","X692.40.small.stuff")



##############
### GATING ###
##############

### NOISE & BEADS
x <- subset(opp, pop==0)
if(draw.gate) plot.cytogram(x, "FSC.small.stuff", "X580.30", main="NOISE & BEADS & SYN")


print("Gating Beads")
if(draw.gate) poly.beads <- getpoly(quiet=TRUE)
beads <- subset(x,inout(x[,c("FSC.small.stuff","X580.30")],poly=poly.beads, bound=TRUE, quiet=TRUE))
opp[row.names(beads),'pop'] <- "beads"


###################
### SAVE PLOT ###
###################
png(paste0(savepath,"/",basename(file),".png"),width=9, height=12, unit='in', res=100)

par(mfrow=c(2,2))
plot.vct.cytogram(opp, "FSC.small.stuff","X692.40.small.stuff")
plot.vct.cytogram(opp, "FSC.small.stuff","X580.30")
plot.vct.cytogram(opp, "X692.40.small.stuff","X580.30")
plot.vct.cytogram(opp, "SSC","X692.40.small.stuff")

dev.off()

###############
### SUMMARY ###
###############

stat.table <- NULL
for(i in unique(opp$pop)){
#print(i)
if(i == 0) next
p <- subset(opp, pop == i)
n <- nrow(p)
if(n ==0) {
fsc <- 0
chl <- 0
    }else{
fsc <- round(median(p$FSC.small.stuff))
chl <- round(median(p$X692.40.small.stuff))
ssc <- round(median(p$SSC))
pe <- round(median(p$X580.30))
}
var <- cbind(i,n,fsc,chl,ssc,pe)
stat.table <- rbind(stat.table, var)
}


table <- data.frame(cbind(stat.table, file=basename(file)))
summary.table <- rbind(summary.table, table)

}

write.csv(summary.table,file=paste(savepath,"/summary.csv", sep=""), row.names=FALSE)



summary <- read.csv(paste0(savepath,"/summary.csv", sep=""))
fr <- read.csv("1.bead_calibration/beadcal_data.csv")
merge <- merge(summary, fr, by.x="file", by.y="File.Name")
merge$abundance <- merge$n/merge$Volume..ul.

par(mfrow=c(1,1), pty='s', cex=1.2)
plot(merge$abundance, merge$beads..ul, log='xy', xlab="Francois", ylab="Sophie", pch=round(merge$inst/100),col=round(merge$Bead.Size..um.*10))
legend('topleft', legend=unique(merge$Bead.Size..um.), pch=1, col=round(merge$Bead.Size..um.*10))
abline(b=1,a=0)

write.csv(merge,file="1.bead_calibration/beadcal_data_FR.csv", row.names=FALSE)
