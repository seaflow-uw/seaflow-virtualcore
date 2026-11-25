#####################################################################
### Bead coordinates (FSC, D1, D2 of optimally positioned beads) ###
### for SeaFlow 130 (second generation, V2)                      ###
#####################################################################
# Same approach as 1.EVT_beadscoordinate.R (740, 751, 989), except that
# gates are drawn automatically (the original used interactive getpoly()),
# and a gating figure is saved for every file in 1.bead_calibration/130_caldata/.
#
#  1. keep particles detected by FSC and aligned on D1 / D2 (|D1 - D2| < width)
#  2. bead cloud in FSC-PE: particles around the bead seed position
#  3. optimally positioned beads = lowest D1 cluster of the bead cloud
#     (beads out of focus sit at higher D1)
#
# Raw data: SeaFlow V2 EVT files from the virtual-core bead experiments
#   ~/Documents/Projects/CALIBRATION/V2_calibration/virtualcore/seaflow/<date>/
# Experiments used: 2021-04-28, 2021-05-26, 2021-06-15, 2021-06-24
# (2021-04-21 excluded: bead D1 far below the other days; 2021-06-09: low quality)
#
# Output: 1.bead_calibration/130-summary.csv (one row per file)
#####################################################################
library(data.table)

repo <- "~/Documents/Codes/seaflow-virtualcore"
path.to.data <- "~/Documents/Projects/CALIBRATION/V2_calibration/virtualcore/seaflow"
setwd(repo)
source(file.path(path.to.data, "functions.R"))   # readSeaflowV2()

experiments <- c("4_28_21", "5_26_21", "6_15_21", "6_24_21")
width <- 2500
dir.create("1.bead_calibration/130_caldata", showWarnings = FALSE)

files <- unlist(lapply(experiments, function(e) list.files(file.path(path.to.data, e), pattern = "micron$", full.names = TRUE)))
# file names use rounded sizes: 035/050/075/1/2/3/6 micron = 0.35/0.53/0.75/1/1.83/3.1/5.7 um beads
SIZES <- c(`035` = 0.35, `050` = 0.53, `075` = 0.75, `1` = 1, `2` = 1.83, `3` = 3.1, `6` = 5.7)
size_of <- function(f) unname(SIZES[sub("micron$", "", sub(".*_", "", basename(f)))])

# Seed positions (median FSC and PE of beads) from the hand-gated bead tables of
# the same experiments (beads_seaflowV2*.csv); the FSC-PE density peak is used
# if a seed is missing.
seed_tab <- rbindlist(lapply(experiments, function(e) {
  f <- list.files(file.path(path.to.data, e), "^beads_seaflowV2b?\\.csv$", full.names = TRUE)
  f <- if (any(grepl("V2b", f))) f[grepl("V2b", f)] else f[1]
  x <- fread(f)[q2.5 == TRUE & pop != "background"]
  x[, size := as.numeric(sub("-.*", "", pop))]
  x[, .(fsc = median(as.numeric(fsc_small)), pe = median(as.numeric(pe))), by = size][, experiment := e]
}))
seed_tab[size == 0.753, size := 0.75][size == 0.99, size := 1]

kde_peak <- function(x, y) {
  k <- MASS::kde2d(x, y, n = 128, lims = c(0, 65535, 0, 65535))
  i <- which(k$z == max(k$z), arr.ind = TRUE)[1, ]
  c(k$x[i[1]], k$y[i[2]])
}

# lowest D1 mode with a valley-to-valley window
low_mode <- function(d1) {
  den <- density(d1, bw = 1200, from = 0, to = 65535, n = 1024)
  pk <- which(diff(sign(diff(den$y))) == -2) + 1
  pk <- pk[den$y[pk] > 0.03 * max(den$y)]
  p <- min(pk)
  vl <- which(diff(sign(diff(den$y))) == 2) + 1
  lo <- if (any(vl < p)) den$x[max(vl[vl < p])] else 0
  hi <- if (any(vl > p)) den$x[min(vl[vl > p])] else 65535
  c(lo, hi)
}

ALL <- NULL
for (f in files) {
  e <- basename(dirname(f)); sz <- size_of(f)
  print(basename(f))
  evt <- readSeaflowV2(f, transform = FALSE)
  evt. <- evt[evt$fsc_small > 1, ]
  aligned <- subset(evt., D2 < D1 + width & D1 < D2 + width)

  sd <- seed_tab[experiment == e & abs(size - sz) < 0.01]
  if (nrow(sd) == 1) { seed <- c(sd$fsc, sd$pe)
  } else {
    cand <- aligned[aligned$fsc_small > 40000 & aligned$pe > 20000, ]
    seed <- kde_peak(cand$fsc_small, cand$pe)
  }
  # bead cloud: ellipse around the seed in FSC-PE
  dfsc <- max(2500, 0.06 * seed[1]); dpe <- max(4000, 0.15 * seed[2])
  b <- aligned[((aligned$fsc_small - seed[1]) / dfsc)^2 + ((aligned$pe - seed[2]) / dpe)^2 <= 1, ]
  b <- subset(b, D1 < 65000 & D2 < 65000)          # remove saturated D1 / D2
  w <- low_mode(b$D1)
  opp <- b[b$D1 >= w[1] & b$D1 <= w[2], ]

  png(paste0("1.bead_calibration/130_caldata/", e, "_", basename(f), ".png"), width = 9, height = 9, units = "in", res = 80)
  par(mfrow = c(2, 2))
  a. <- if (nrow(aligned) > 1e5) aligned[sample(nrow(aligned), 1e5), ] else aligned
  smoothScatter(a.$fsc_small, a.$pe, xlim = c(0, 65535), ylim = c(0, 65535), nrpoints = 0,
                xlab = "fsc_small", ylab = "pe", main = paste(e, basename(f)))
  th <- seq(0, 2 * pi, length.out = 100)
  lines(seed[1] + dfsc * cos(th), seed[2] + dpe * sin(th), col = 2, lwd = 2)
  smoothScatter(b$fsc_small, b$D1, xlim = c(0, 65535), ylim = c(0, 65535), nrpoints = 0, xlab = "fsc_small", ylab = "D1",
                main = "bead cloud; red = optimally positioned")
  abline(h = w, col = 2, lwd = 2)
  smoothScatter(b$fsc_small, b$D2, xlim = c(0, 65535), ylim = c(0, 65535), nrpoints = 0, xlab = "fsc_small", ylab = "D2")
  points(opp$fsc_small, opp$D2, pch = ".", col = 2)
  hist(b$D1, breaks = 100, main = "D1 of bead cloud", xlab = "D1"); abline(v = w, col = 2, lwd = 2)
  dev.off()

  q <- function(x, p) as.numeric(quantile(x, p))
  ALL <- rbind(ALL, data.frame(file = paste0("130_caldata/", e, "_", basename(f)), experiment = e, size = sz,
                               beads = nrow(opp),
                               fsc.med = q(opp$fsc_small, .5), fsc.low = q(opp$fsc_small, .025), fsc.high = q(opp$fsc_small, .975),
                               d1.med = q(opp$D1, .5), d1.low = q(opp$D1, .025), d1.high = q(opp$D1, .975),
                               d2.med = q(opp$D2, .5), d2.high = q(opp$D2, .975), d2.low = q(opp$D2, .025)))
}

write.csv(ALL, "1.bead_calibration/130-summary.csv", quote = FALSE, row.names = FALSE)
