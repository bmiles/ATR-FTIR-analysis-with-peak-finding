library('ProjectTemplate')
load.project()

#commands to get peaks to work
library("Peaks")
library.dynam('Peaks', 'Peaks', lib.loc=NULL)


source('src/func.R')
groupHooks <- c('sample3', 'sample4') #sample is the optics bacgrkound measurement
groupColors <- brewer.pal(length(groupHooks), "Dark2")

setwd('data')
files2SpectraObject(gr.crit = groupHooks, gr.cols = groupColors,
                    freq.unit = "wavenumber cm^-1",
                    int.unit = "Absorbance Units",
                    descrip = "ATR-FTIR data",
                    format = "csv",
                    alignTMS = FALSE,
                    out.file = "specData",
                    debug = FALSE)
setwd('..')
load("data/specData.RData")
spec <- saveLoadReference
baseSpec <- baselineSpec(saveLoadReference, int=FALSE, retC=TRUE, method="rfbaseline")

# remove noisy, uninteresting region:
newIR <- removeFreq(baseSpec, rem.freq =
                      baseSpec$freq > 1800 & baseSpec$freq < 2500)

## Set ATR crystal Material. 'Ge' or 'Diamond'
crystalRange <- getCrystalRangeFor('Ge')
# ## Print Backgrounds pre and post purge
# newFileName <- paste(getwd(), "/graphs/Background spectra pre and post purge with nitrogen, 17-06-2014", ".pdf", sep = "")
#   pdf(file=newFileName, paper="a4r", width=9)
#   par(mfrow = c(1,1), mar = c(5.1,5.1,4.1,2.1))
#   plotSpectra(baseSpec,
#               which=c(1,3),
#               xlim=crystalRange, #Ge ATR crystal range.
#               title="Background spectra pre and post purge with nitrogen, 17-06-2014",
#               offset= 1,
#               yrange=c(-10,7),
#               xaxt = "n")
#   at <- seq(from = 0, to = max(baseSpec$freq), by = 100)
#   axis(side = 1, at = at, las = 2, hadj = 0.9)
# dev.off()

newFileName <- paste(getwd(), "/graphs/ATR-FTIR mineralized samples, 14-06-2014", ".pdf", sep = "")
  pdf(file=newFileName, paper="a4r", width=9)
  par(mfrow = c(1,1), mar = c(5.1,5.1,4.1,2.1))
  plotSpectra(baseSpec,
              which=c(6:1),
              xlim=crystalRange, #Ge ATR crystal range.
              title="ATR-FTIR mineralized samples, 14-06-2014",
              offset= 0.01,
              yrange=c(-0.01,0.25),
              xaxt = "n")
  at <- seq(from = 0, to = max(baseSpec$freq), by = 100)
  axis(side = 1, at = at, las = 2, hadj = 0.9)
dev.off()

newFileName <- paste(getwd(), "/graphs/ATR-FTIR mineralized samples, 14-06-2014, detail", ".pdf", sep = "")
  pdf(file=newFileName, paper="a4r", width=9)
  par(mfrow = c(1,1), mar = c(5.1,5.1,4.1,2.1))
  plotSpectra(baseSpec,
              which=c(6:1),
              xlim=c(1800,800),
              title="ATR-FTIR mineralized samples, 14-06-2014, detail",
              offset= 0.01,
              yrange=c(0,0.25),
              xaxt = "n",
              lab.pos = 1580)
  at <- seq(from = 0, to = max(baseSpec$freq), by = 20)
  axis(side = 1, at = at, las = 2, hadj = 0.9)
dev.off()

peakFinder(2)
ggsave('graphs/peaks for spectra 2.pdf', width=14)
