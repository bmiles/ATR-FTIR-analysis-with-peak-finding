gplotSpectra <- function(spectra, which = c(2),
  yrange = range(spectra$data),
  offset = 0.0, amplify = 1.0,
  lab.pos = mean(spectra$freq), ...) {

# Function to plot multiple spectra @ specified expansions & decorate
# Part of the ChemoSpec package
# Bryan Hanson, DePauw University, June 2009

  if (missing(spectra)) stop("No spectral data provided")
  chkSpectra(spectra)

  # set up and plot the first spectrum

  spectrum <- spectra$data[which[1],]*amplify
  yseq <- seq(from= yrange[1], to=yrange[2])
  #ggplot
   plots <- ggplot() + geom_path(aes(x= spectra$freq, y=spectrum))
     #
    #  which <- which[-1] # first spectrum already plotted so remove it from the list
    #  count <- 0 # get the other spectra and plot them as well
    #  for(n in which) {
    #    count <- count + 1
    #    spectrum <- (spectra$data[n,]+(offset*count))*amplify
    #    #points(spectra$freq, spectrum, type = "l", col = spectra$colors[n], ...)
    #    plots <- plots + geom_path(aes(x= spectra$freq, y=spectrum, colour= spectra$colors[n]))
    #    #lab.y <- spectrum[spec.index]
    #    #text(lab.x, lab.y, labels = spectra$names[n], pos = 3, cex = 0.75)
    #    }
     #geom_path(aes(x= spec$freq, y=spec$data[3,]))

  #qplot

  # qplot(spectra$freq, spectrum, geom = "path",
  # xlab = spectra$unit[1], ylab = spectra$unit[2],
  # ylim = yrange, col = spectra$colors[which[1]],
  # frame.plot = FALSE, ...)
  # grid(ny = NA, lty = 1)
  # lab.x <- lab.pos
  # spec.index <- findInterval(lab.x, sort(spectra$freq))
  # lab.y <- spectrum[spec.index]
  # text(lab.x, lab.y, labels = spectra$names[which[1]], pos = 3, cex = 0.75)
  #
  # which <- which[-1] # first spectrum already plotted so remove it from the list
  # count <- 0 # get the other spectra and plot them as well
  # for(n in which) {
  # 	count <- count + 1
  # 	spectrum <- (spectra$data[n,]+(offset*count))*amplify
  # 	points(spectra$freq, spectrum, type = "l", col = spectra$colors[n], ...)
  # 	lab.y <- spectrum[spec.index]
  # 	text(lab.x, lab.y, labels = spectra$names[n], pos = 3, cex = 0.75)
  # 	}
}
# Heat Map, returns a heat map of intensities, takes a wavenumber or range (and calculates the mean)
irHeatmap <- function (wavenumber) {
  # Setup these parametric variables.
  pH <- c(3,7.3,12)
  collagenConc <- c(0,10,30)

  # Array to hold the intensities for a given wavenumber
  absVal <- c()


  # Loop gets intensities.
  for (i in c(1:length(newIR$groups))) {
    nextVal <-mean(c(newIR$data[i, wavenumber]))
    absVal <- c(absVal, nextVal)
  }

  # creates a matrix of the intensity values.
  vnames <-c()
  for (i in c(1:length(newIR$names))) {
   nextName <-substring(c(newIR$names), 1, 2)
   vnames <- c(vnames, nextName)
  }
  mnames <- matrix(vnames, nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(pH, collagenConc))

  mdat <- matrix(c(absVal), nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(pH, collagenConc))
  colors <- brewer.pal(9, "YlGnBu")
  my_palette <- colorRampPalette(colors)(n = 299)
  heatmap.2(mdat,
            main = "Collagen deposition dependence on pH and concentration", # heat map title
            cellnote = mnames,
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(12,12),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="none",     # only draw a row dendrogram
            Colv="NA",            # turn off column clustering
            Rowv="NA",
            xlab = "Collagen Concentration v/v %",
            ylab = "pH",
            srtCol = 0,
            offsetRow = 2,
            offsetCol = 3)
}
getCrystalRangeFor <- function (crystal) {
  if (crystal == 'Ge') {
    c(4400,780)
  }
  else if (crystal == 'Diamond') {
    c(4400,400)
  }
  else print(crystal)
}
peakFinder <- function(specVector) {
  #specVectorExample <- numeric 1 or 3 etc.
  freq <- baseSpec$freq
  intData <- baseSpec$data[specVector,]
  decon <- SpectrumSearch(intData*1000,sigma=3.0,threshold=7, background=FALSE, iterations=20, markov=FALSE,window=3)

  peaksLoc <- decon$pos
  peaks <- freq[peaksLoc]


  df <- data.frame(freq, intData)
  p <- ggplot(df, aes(x=freq, y=intData))+
    geom_path() +
    scale_x_reverse(breaks = signif(seq(min(freq), max(freq), by = 200)))
  # Note the xend value being multiplied by -1, this was required when the scale was reversed.
  a <- annotate("segment",y=min(intData),yend=max(intData), x=peaks, xend= peaks * -1, color="red")
  t <- annotate("text", x = peaks, y =max(intData) * 1.05, label = peaks, angle=90, color="red")
  glabs <- labs(title="FTIR spectra with peaks indicated",
              x = "wavenumber cm^-1",
              y = "Absorbance Units")
  p + a + t + glabs

}
irHeatmapRMS <- function (wavenumber) {
  # Setup these parametric variables.
  pH <- c(3,7.3,12)
  collagenConc <- c(0,10,30)

  # Array to hold the intensities for a given wavenumber
  absVal <- c()


  # Loop gets intensities.
  for (i in c(1:length(newIR$groups))) {
    y <- c(newIR$data[i, wavenumber])
    nextVal <-sqrt(sum(y^2)/length(y))
    absVal <- c(absVal, nextVal)
  }

  # creates a matrix of the intensity values.
  vnames <-c()
  for (i in c(1:length(newIR$names))) {
    nextName <-substring(c(newIR$names), 1, 2)
    vnames <- c(vnames, nextName)
  }
  mnames <- matrix(vnames, nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(pH, collagenConc))

  mdat <- matrix(c(absVal), nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(pH, collagenConc))
  colors <- brewer.pal(9, "YlGnBu")
  my_palette <- colorRampPalette(colors)(n = 299)
  heatmap.2(mdat,
            main = "Collagen deposition dependence on pH and concentration", # heat map title
            cellnote = mnames,
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(12,12),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="none",     # only draw a row dendrogram
            Colv="NA",            # turn off column clustering
            Rowv="NA",
            xlab = "Collagen Concentration v/v %",
            ylab = "pH",
            srtCol = 0,
            offsetRow = 2,
            offsetCol = 3)
}
