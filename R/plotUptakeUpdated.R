plotUptakeUpdated <- function(pep1,pep2,times,pepIDs = peptide.identifications,pepnum){

    # Modify tables
    new1 <- pep1
    new2 <- pep2


    maxYvalue <- pepIDs[pepnum,15]
    maxYlim <- max(max(max(new1),max(new2)),maxYvalue)

    times[1] <- 1
    plot(times,new1,type='b',xlab='Time (min)',log = 'x',ylab = 'Relative Deuterium Uptake (Da)', main = paste("Deuterium Level Peptide Number",as.character(pepnum),sep=" "),col="red",pch=4,ylim=c(0,maxYvalue+10))
    lines(times,new2,type='b',col="green",pch=3)
    abline(h = maxYvalue,col="blue",lty =2)
    text(4, maxYvalue +3, "Theoretical Maximum Uptake", col = "blue",cex = 0.75)

    legend('topleft',c('unbound','bound'),
           lty=1,lwd=c(1,1),col=c('red','green'))

}
