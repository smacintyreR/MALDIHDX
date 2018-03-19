PlotUptakeCompare <- function(pepnum,pepIDS = peptide.identifications,AllCents,times,Conditions = c("Bound","Unbound")){

    
  NoRep <- length(AllCents)/2
  CentsCondit1 <- AllCents[1:NoRep]
  CentsCondit2 <- AllCents[(NoRep+1):2*NoRep]
  
  AvCentsCondit1 <- AvCent(CentsCondit1)
  AvCentsCondit2 <- AvCent(CentsCondit2)

  
  times <- as.numeric(as.vector(CentsCondit1[[1]][,1]))
  times[1] <- 1
  noTimes <- length(times)


  # Add the theoretical undeuterated
  #UnCent <- pepIDS[pepnum,18]

  # HDX Uptake values
  #CentsCond1Mod <- rbind(as.data.frame(CentsCondit1[,2]),UnCent)

  UptakesCond1 <- AvCentsCondit1 - AvCentsCondit1[1,]

  # HDX Uptake values
  
  
  #CentsCond2Mod <- rbind(as.data.frame(CentsCondit2[,2]),UnCent)

  UptakesCond2 <- AvCentsCondit2 - AvCentsCondit2[1,]

  xmax <- times[noTimes] + 3
  YmaxTh <- pepIDS[pepnum,15]

  YvalueMax <- max(YmaxTh,max(UptakesCond2))
  # Plotting
  plot(times,UptakesCond1[,1],log='x',ylim = c(0,YvalueMax),col = "blue",pch=4,ylab="Relative Deuterium Uptake (AMU)",xlab = "Time (seconds)",main = paste("Back Exchange Peptide",as.character(pepnum),"Flag Binding",sep="_"),type="b")
  points(times,UptakesCond2[,1],col="red",pch=3,type="b")

  abline(h = YmaxTh,col="blue",lty =2)
  text(times[noTimes]-3555, YmaxTh - 0.25, "Theoretical Maximum Uptake", col = "blue",cex = 0.75)

  legend('topright',Conditions,
         pch=c(3,4),col=c('red','blue'))

}
