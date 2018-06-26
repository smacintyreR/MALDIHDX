#' Average Centroid Calculation
#' @description  Averages the centroids for each time point across each
#' technical replicate
#'
#' @param AllpeptideSamples A singular mass spectrum
#' @keywords centroids,average,sample
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' [a,b,envelope] <- widthFinder(peaks,spectrum,bpint = 0.5)


AvCent <- function(AllpeptideSamples){
    
    numSam <- length(AllpeptideSamples)
    numTime <- dim(AllpeptideSamples[[1]])[1]
    
    AvVec <- rowMeans(cbind(AllpeptideSamples[[1]][,2],AllpeptideSamples[[2]][,2],AllpeptideSamples[[3]][,2]))
    AvVec <- as.data.frame(AvVec)
    
    
    colnames(AvVec)<-"Average Centroid"
    
    
    return(AvVec)
    
    
}





#' Manual Peak Selection
#' @description  Enables the locator to be used to select peaks and hence
#' calculates the centroid and peak width
#'
#'
#' @param spectrum A singular mass spectrum
#' @param i sample number
#' @param j time point
#' @param TimeP timepoints
#' @param mass masses
#' @param bpint base peak intensity
#' @keywords peaks, manual, locator
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' [c,v] <- manualCentroids2(s,1,2,mass)


manualCentroids2 <- function(spectrum,i,j,TimeP,mass,bpint=0.5){
    
    plot(spectrum,main=paste("Sample",j,"Time",TimeP[i],sep=" "))
    p <- locator(type="p",col="red",pch=1)
    p <- as.data.frame(p)
    p <- p[order(p$x),]
    p <- createMassPeaks(mass=p$x,intensity = p$y)
    #cat("\n","Enter peak width (as % of base peak intensity):","\n")
    #bpint <- scan(n=1)
    v <- widthFinder(p,spectrum,bpint)
    w <- v[2]-v[1]
    c <- centroidCalc(v,spectrum)
    
    interp <- linearInterp(p)
    
    plotLin(interp,p)
    plotWidth(v)
    centroidPlot(c,spectrum)
    legend('topright',c('Spectrum','Distribution Width','Centroid','Peak envelope'),
           lty=1,lwd=c(1,4,4,1),col=c('black','orange','green','red'))
    
    cat ("Press [enter] to continue")
    line <- readline()
    
    
    
    return(c(c,w))
    
}



#' Centroid relative uptakes calculation and export for a single peptide
#' @description  Calculates the centroids for each sample and timepoint for
#' each
#'
#'
#' @param spectrum A singular mass spectrum
#' @param i sample number
#' @param j time point
#' @param TimeP timepoints
#' @param mass masses
#' @param bpint base peak intensity
#' @keywords peaks, manual, locator
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' [c,v] <- manualCentroids2(s,1,2,mass)






mainCentNewMod <- function(samplematrix,charge=1,SNR=10,bpint=0.5,sample=0,time=0){
    
    Samples <- unique(samplematrix[,1])
    NoSamples <- length(Samples)  #Equal to 3 here
    w <- c() # vector to hold widths
    cp <- c() #Vector to hold positional centroid
    c <- c()
    TimeP <- unique(samplematrix[,2]) # The time points for this experiment
    noTimeP <- length(TimeP) #Number of time points
    new <- samplematrix[,-c(1,2)]
    masses <- as.numeric(colnames(new))
    
    
    out <- list()
    curDir <- getwd()
    dir.create(paste("Peptide",masses[1],"Centroid","Plots",sep="_"))
    setwd(paste("Peptide",masses[1],"Centroid","Plots",sep="_"))
    
    
    
    for (j in 1:NoSamples){
        
        for(i in 1:noTimeP){
            
            index <- i+(j-1)*noTimeP
            s <- createMassSpectrum(mass=masses,intensity=as.numeric(new[index,]))
            p <- peakPick(s,SNR)
            v <- widthFinder(p,s,bpint)
            w[i] <- (v[2]-v[1])*charge
            cp[i] <- centroidCalc(v,s)
            
            c[i] <- cp[i]
            
            interp <- linearInterp(p)
            
            pdf(paste(masses[1],"CentroidPlot","Sample",j,"-Time",TimeP[i],".pdf",sep=''))
            plot(s,main=paste("Sample",j,"Time",TimeP[i],sep=" "))
            plotLin(interp,p)
            plotWidth(v)
            centroidPlot(cp[i],s)
            legend('topright',c('Spectrum','Distribution Width','Centroid','Peak envelope'),
                   lty=1,lwd=c(1,4,4,1),col=c('black','orange','green','red'))
            dev.off()
            
        }
        summary <- cbind(TimeP,c)
        summary <- as.data.frame(summary)
        summary[,2] <- as.numeric(as.character(summary[,2]))
        summary[,3] <- summary[,2] - summary[1,2]
        summary[,4] <- w
        colnames(summary) <- c("time (min)","centroid (Da)","Rel D Lvl (Da)","width (Da)")
        out[[j]] <- summary
    }
    
    setwd(curDir)
    
    dir.create(paste("Peptide",masses[1],"Uptake","Plots",sep="_"))
    setwd(paste("Peptide",masses[1],"Uptake","Plots",sep="_"))
    
    noRep <- NoSamples/2  # Assuming a two state experiment
    
    Unbound <- out[1:noRep]
    Bound <- out[(noRep+1):NoSamples]
    
    times <- out[[1]][,1]
    
    AvCentsUnbound <- AvCent(Unbound)
    AvCentsBound <- AvCent(Bound)
    
    png(paste(masses[1],"UptakePlot1",".png",sep=''))
    
    plotUptakeUpdated(AvCentsUnbound[[1]],AvCentsBound[[1]],times)
    
    dev.off()
    
    setwd(curDir)
    
    
    
    
    return(out)
    
}




#' Centroid calculation and uptakes for a complete peptide
#' @description  Allows the calculation of centroids for each time point and
#' sample for a single peptide. Also prompts manual centroid calculation if
#' automatic value is not appropriate. Each centroid is saved as a pdf in the
#' current directory.
#'
#'
#'
#' @param samplematrix Matrix of all timepoints and samples for a single peptide
#' @param charge
#' @param SNR time point
#' @param bpint timepoints
#' @param sample masses
#' @param time base peak intensity
#' @keywords peaks, manual, locator
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' Centroids <- mainCentMerge(Peptide1,SNR=5)


mainCentMerge <- function(samplematrix,charge=1,SNR=10,bpint=0.5,sample=0,time=0){
    
    Samples <- samplematrix[,1]
    SamplesU <- unique(Samples) #List of Samples (E.g A1,A2...)
    NoSamples <- length(SamplesU)  #Number of technical replicates
    w <- c() # vector to hold widths
    cp <- c() #Vector to hold positional centroid
    c <- c() # Vector to h0ld reported centroid
    flag <- c() # vector to hold flagged centroid data
    TimeP <- unique(samplematrix[,2]) # The time points for this experiment
    noTimeP <- length(TimeP) # Number of time points
    new <- samplematrix[,-c(1,2)] # Remove column 1&2
    masses <- as.numeric(colnames(new)) # Extract masses
    
    # Initialise list for summary tables for each replicate
    out <- list()
    
    
    
    for (j in 1:NoSamples){
        
        for(i in 1:noTimeP){
            
            index <- i+(j-1)*noTimeP
            s <- createMassSpectrum(mass=masses,intensity=as.numeric(new[index,]))
            p <- peakPick(s,SNR=SNR)
            
            if(length(p) > 1){
                v <- widthFinder(p,s,bpint)
                w[i] <- (v[2]-v[1])*charge
                cp[i] <- centroidCalc(v,s)
                
                c[i] <- cp[i]
                
                interp <- linearInterp(p)
                
                #png(filename = paste(masses[1],"CentroidPlot","Sample",j,"-Time",TimeP[i],".png",sep=''))
                plot(s,main=paste("Sample",Samples[index],"Time",TimeP[i],sep=" "))
                plotLin(interp,p)
                plotWidth(v)
                centroidPlot(cp[i],s)
                legend('topright',c('Spectrum','Distribution Width','Centroid','Peak envelope'),
                       lty=1,lwd=c(1,4,4,1),col=c('black','orange','green','red'))
                cat ("Press [Y] if satisfied with centroid, [N] if require manual editing, [F] to flag centroid: ")
                line <- readline()
                if(line=="Y"){
                    flag[i] <- F
                }
                
                else if(line == "F"){
                    flag[i] <- T  
                }
                
                
                else{
                    c[i] <-manualCentroids2(s,i,j,TimeP,masses)[1]
                    flag[i] <- F
                }
            }
            
            else{
                plot(s,main=paste("Sample",Samples[index],"Time",TimeP[i],sep=" "))
                cat ("This centroid will be flagged, press [enter] to continue:")
                line <- readline()
                flag[i] <- T
                c[i] <- 0
                w[i] <- 0
            }
            #dev.off()
            
        }
        summary <- cbind(TimeP,c)
        summary <- as.data.frame(summary)
        summary[,2] <- as.numeric(as.character(summary[,2]))
        summary[,3] <- summary[,2] - summary[1,2]
        summary[,4] <- w
        summary[,5] <- flag
        colnames(summary) <- c("time (min)","centroid (Da)","Rel D Lvl (Da)","width (Da)","FlAGGED")
        out[[j]] <- summary
    }
    
    
    return(out)
    
}







mainCentNewMod2 <- function(samplematrix,charge=1,SNR=5,bpint=0.5){
    
    Samples <- unique(samplematrix[,1])
    NoSamples <- length(Samples)  #Equal to 3 here
    w <- c() # vector to hold widths
    cp <- c() #Vector to hold positional centroid
    c <- c()
    TimeP <- unique(samplematrix[,2]) # The time points for this experiment
    noTimeP <- length(TimeP) #Number of time points
    new <- samplematrix[,-c(1,2)]
    masses <- as.numeric(colnames(new))
    
    
    out <- list()
    
    for (j in 1:NoSamples){
        
        for(i in 1:noTimeP){
            
            index <- i+(j-1)*noTimeP
            s <- createMassSpectrum(mass=masses,intensity=as.numeric(new[index,]))
            p <- peakPick(s,SNR)
            
            if(length(p)>1){
                v <- widthFinder(p,s,bpint)
                w[i] <- (v[2]-v[1])*charge
                c[i] <- centroidCalc(v,s)
            }
            
            else{
                w[i] <- NA
                c[i] <- NA
            }
            
            
            
            
        }
        
        summary <- cbind(TimeP,c)
        summary <- as.data.frame(summary)
        summary[,2] <- as.numeric(as.character(summary[,2]))
        summary[,3] <- summary[,2] - summary[1,2]
        summary[,4] <- w
        colnames(summary) <- c("time (min)","centroid (Da)","Rel D Lvl (Da)","width (Da)")
        out[[j]] <- summary
    }
    
    return(out)
    
}
