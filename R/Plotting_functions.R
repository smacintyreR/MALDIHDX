# Plotting functions for HX - Macspress

#' Linear Interpolation Plotting Function
#' @description  Given the linear interpolation data and the mass peaks,
#' plots both over the spectrum
#' @param peaks list of mass peaks
#' @keywords linear,interpolation, plot, peaks
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' plotLin(interp,peaks)
plotLin <- function(interp,peaks){
  points(peaks,col="red",lwd=4)
  lines(interp,col="red")
}

# Width segment

plotWidth <- function(v){
  a <- v[1]
  b <- v[2]
  envelope <- v[3]

  segments(a,envelope,b,envelope,col="orange",lwd=4)
}

# Centroid segment

centroidPlot <- function(cent,s){
  segments(cent,0,cent,max(intensity(s)),col="green",lwd=4)
}

# Legend



# Plots the Deuterium uptake curves for the two states
plotUptake <- function(pep1,pep2){

  # Modify tables
  new1 <- pep1
  new2 <- pep2
  new1[,1] <- as.numeric(as.character(new1[,1]))
  new2[,1] <- as.numeric(as.character(new2[,1]))
  new1 <- new1[-1,]
  new2 <- new2[-1,]
  new1 <- new1[-dim(new1)[1],]
  new2 <- new2[-dim(new2)[1],]

  plot(new1[,1],new1[,3],type='b',xlab='Time (min)',ylab = 'Relative Deuterium Uptake (Da)', main = "Deuterium Level",log="x",col="red",pch=4,ylim=c(0,14))
  lines(new2[,1],new2[,3],type='b',col="green",pch=3)

  legend('topleft',c('mutant','wild type'),
         lty=1,lwd=c(1,1),col=c('red','green'))

}


plotUptake <- function(pep1,pep2){

  # Modify tables
  new1 <- pep1
  new2 <- pep2
  new1[,1] <- as.numeric(as.character(new1[,1]))
  new2[,1] <- as.numeric(as.character(new2[,1]))
  new1 <- new1[-1,]
  new2 <- new2[-1,]
  new1 <- new1[-dim(new1)[1],]
  new2 <- new2[-dim(new2)[1],]

  plot(new1[,1],new1[,3],type='b',xlab='Time (min)',ylab = 'Relative Deuterium Uptake (Da)', main = "Deuterium Level",log="x",col="red",pch=4,ylim=c(0,14))
  lines(new2[,1],new2[,3],type='b',col="green",pch=3)

  legend('topleft',c('mutant','wild type'),
         lty=1,lwd=c(1,1),col=c('red','green'))

}

# Uses all of the plot functions to plot and save all of the centroided, peakwidthed spectra
# Inefficient at the momemt given it recalculates centroids etc
# Will develop further
# Does one peptide at a time
mainCentPlots <- function(file){



  L <- dim(file)[2]/2

  file <- file[complete.cases(file),]

  for(i in 1:L){

    s <- createMassSpectrum(mass=file[,2*i-1],intensity=file[,2*i])
    p <- peakPick(s)
    v <- widthFinder(p,s)
    c <- centroidCalc(v,s)
    interp <- linearInterp(p)

    png(filename = paste("CentroidPlot",i,".png",sep=''))
    plot(s)
    plotLin(interp,p)
    plotWidth(v)
    centroidPlot(c,s)
    dev.off()
  }

}
