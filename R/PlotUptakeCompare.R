PlotUptakeCompare <- function(pepnum,pepIDS = peptide.identifications,CentsCondit1,CentsCondit2,times){


  noTimes <- length(times)


  # Add the theoretical undeuterated
  UnCent <- pepIDS[pepnum,18]

  # HDX Uptake values
  CentsCond1Mod <- rbind(as.data.frame(CentsCondit1[,2]),UnCent)

  UptakesCond1 <- CentsCond1Mod - CentsCond1Mod[7,]

  # HDX Uptake values
  CentsCond2Mod <- rbind(as.data.frame(CentsCondit2[,2]),UnCent)

  UptakesCond2 <- CentsCond2Mod - CentsCond2Mod[7,]

  xmax <- times[noTimes] + 3
  YmaxTh <- pepIDS[pepnum,15]

  YvalueMax <- max(YmaxTh,max(UptakesCond2))
  # Plotting
  plot(times,UptakesCond1[,1],log='x',ylim = c(0,YvalueMax),col = "blue",pch=4,ylab="Relative Deuterium Uptake (AMU)",xlab = "Time (minutes)",main = paste("Back Exchange Peptide",as.character(pepnum),"Temperature Dependence",sep="_"))
  points(times,UptakesCond2[,1],col="red",pch=3)

  abline(h = YmaxTh,col="blue",lty =2)
  text(times[noTimes]-60, YmaxTh - 0.25, "Theoretical Maximum Uptake", col = "blue",cex = 0.75)

  legend('bottomleft',c(expression(paste("0 ",degree,"C")),expression(paste("- 20 ",degree,"C"))),
         pch=c(4,3),col=c('red','blue'))




}
