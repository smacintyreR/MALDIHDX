PlotUptakeCompare <- function(pep.num, pep.ids = peptide.identifications,
                              all.cents, times,
                              conditions = c("Bound","Unbound")){

  # Separate centroids based on experimental conditions
  rep.no <- length(all.cents)/2
  cents.condit.1 <- all.cents[1:rep.no]
  cents.condit.2 <- all.cents[ (rep.no + 1):2*rep.no]
    
  # Average centroids across technical replicates
  av.cents.condit.1 <- AvCent(cents.condit.1)
  av.cents.condit.2 <- AvCent(cents.condit.2)
    
  # Extract time points
  times <- as.numeric(as.vector(cents.condit.1[[1]][, 1]))
  times[1] <- 1
  no.times <- length(times)
    
  # Calculate uptake values from centroids
  uptakes.cond.1 <- av.cents.condit.1 - av.cents.condit.1[1, ]
  uptakes.cond.2 <- av.cents.condit.2 - av.cents.condit.2[1, ]
       
  # Set plot dimensions
  x.max <- times[no.times] + 3
  y.max.th <- pep.ids[pep.num, 15]
  y.value.max <- max(y.max.th, max(uptakes.cond.2))
    
  # Plot the uptake curves 
  plot(times, uptakes.cond.1[, 1], log = 'x',ylim = c(0, y.value.max), 
       col = "blue", pch = 4, ylab = "Relative Deuterium Uptake (AMU)",
       xlab = "Time (seconds)", main = paste("Back Exchange Peptide",
    as.character(pepnum), "Flag Binding", sep="_"), type = "b")
  points(times,uptakes.cond.2[,1],col="red",pch=3,type="b")
  abline(h = y.max.th, col="blue", lty=2)
    
  # Add label and legend
  text(times[no.times] - 3555, y.max.th - 0.25, "Theoretical Maximum Uptake", 
  col = "blue", cex = 0.75)
  legend('topright', conditions, pch=c(3,4), col=c('red','blue'))
}
