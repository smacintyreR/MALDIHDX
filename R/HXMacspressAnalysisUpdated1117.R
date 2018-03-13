# Analysis functions for HX-Macspress


#' Peak Picking Function
#'
#' Identifies peaks in the spectrum
#' @param spectrum A singular mass spectrum
#' @keywords peaks, SNR
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' p <- peakPick(s)
peakPick <- function(spectrum,SNR=10){


  peaks <- detectPeaks(spectrum,SNR=SNR)

  return(peaks)

}


# Finds the peak width using two different types of linear interpolation
# Interpolates between the peaks if the peak envelope is below the threshold
# Default base peak intensity % = 20

widthFinder <- function(peaks,spectrum,bpint =0.2){

  # Calculate the intensity cut-off
  envelope = bpint*max(intensity(spectrum))

  #Find which (x1,y1) and (x2,y2) to use

  L <- length(mass(peaks))
  M <- length(mass(spectrum))
  # Set (x1,y1) and (x2,y2) to peak 1 and 2

  x1 <- mass(peaks)[1]
  y1 <- intensity(peaks)[1]
  x2 <- mass(peaks)[2]
  y2 <- intensity(peaks)[2]

  # Case 1: first peak lies above threshold
  # Must interpolate using spectrum
  if(y1 > envelope){

	# Set initial points as the first/second in spectrum
    x1 <- mass(spectrum)[1]
    y1 <- intensity(spectrum)[1]
    x2 <- mass(spectrum)[2]
    y2 <- intensity(spectrum)[2]

	# Search through pairs of points until the envelope
	# value is bounded between y1 ad y2
    for(i in 1:M){
      if(y2 < envelope){
        x2 <- mass(spectrum)[i+2]
        y2 <- intensity(spectrum)[i+2]
        x1 <- mass(spectrum)[i+1]
        y1 <- intensity(spectrum)[i+1]
      }

      else{

		# Now calculate the left window by interpolation
        a <- linearsolve(c(x1 ,y1),c(x2, y2),envelope)
      }
    }

  }

  # Case 2: first peak lies below threshold
  # Interpolate using peaks
  else{
    for(i in 1:L){
      if(y2 < envelope){
        x2 <- mass(peaks)[i+2]
        y2 <- intensity(peaks)[i+2]
        x1 <- mass(peaks)[i+1]
        y1 <- intensity(peaks)[i+1]
      }

      else{

        a <- linearsolve(c(x1 ,y1),c(x2, y2),envelope)
      }
    }
  }


  # Finding the right end point with a symmetrical process
  x1 <- mass(peaks)[L]
  y1 <- intensity(peaks)[L]
  x2 <- mass(peaks)[L-1]
  y2 <- intensity(peaks)[L-1]

  if(y1 > envelope){

    x1 <- mass(spectrum)[L]
    y1 <- intensity(spectrum)[L]
    x2 <- mass(spectrum)[L-1]
    y2 <- intensity(spectrum)[L-1]

    for(i in M:1){
      if(y2 < envelope){
        x2 <- mass(spectrum)[i-2]
        y2 <- intensity(spectrum)[i-2]
        x1 <- mass(spectrum)[i-1]
        y1 <- intensity(spectrum)[i-1]
      }

      else{
        b <- linearsolve(c(x1, y1),c(x2, y2),envelope)
      }
    }

  }

  else{

  for(i in L:1){
    if(y2 < envelope){
      x2 <- mass(peaks)[i-2]
      y2 <- intensity(peaks)[i-2]
      x1 <- mass(peaks)[i-1]
      y1 <- intensity(peaks)[i-1]
    }

    else{
      b <- linearsolve(c(x1, y1),c(x2, y2),envelope)
    }
  }
  }

  return(c(a,b,envelope))

}





centroidCalc <- function(v,spectrum){

  a <- v[1]
  b <- v[2]

  new <- subset(spectrum,mass(spectrum) < b & mass(spectrum) > a)

  cent <- sum(intensity(new)*mass(new))/sum(intensity(new))

  return(cent)

}






# Main function (calculates all centroids, Relative uptake and width)

# Note: Takes HX-Express type input file for a time course for
# one peptide, a list of the time points and the charge of the peptide

mainCent <- function(file,times,charge){


  # Calculate the number of peptides
  L <- dim(file)[2]/2

  # Empty vector of peak widths
  w <- c()

  # Removes rows which have missing data
  file <- file[complete.cases(file),]

  for(i in 1:L){

    s <- createMassSpectrum(mass=file[,2*i-1],intensity=file[,2*i])
    p <- peakPick(s)
    v <- widthFinder(p,s)
    w[i] <- (v[2]-v[1])*charge
	# Where does 1.00794 come from??
    c[i] <- centroidCalc(v,s)

  }

  # Output in form of HX Express output
  new <- cbind(times,c)
  new <- as.data.frame(new)
  new[,2] <- as.numeric(as.character(new[,2]))
  new[,3] <- new[,2] - new[1,2]
  new[,4] <- w
  colnames(new) <- c("time (min)","centroid (Da)","Rel D Lvl (Da)","width (Da)")
  return(new)

}



# Function solves a linear equation at the
# value = value and passing through the points
# p1 and p2
linearsolve <- function(p1,p2,value){
  x1 <- p1[1]
  y1 <- p1[2]
  x2 <- p2[1]
  y2 <- p2[2]

  m <- (y2-y1)/(x2-x1)

  x <- (value - y1 +m*x1)/m

  return(x)
}



# Function calculates the linear interpolation between
# the peaks, needed for plotting tools
linearInterp <- function(peaks){

  peaksmod <-cbind(mass(peaks),intensity(peaks))

  interp <- approx(x=peaksmod[,1],y=peaksmod[,2],method="linear",n=100,rule=1)
  return(interp)

}





# mainCentMerge(samplematrix,charge=1,SNR=10,bpint=0.5,sample=0,time=0)
# Author: Sam MacIntyre
# Input: The feature matrix for a single peptide in the form of a large matrix
# containing all of the technical replicates
# Output: Outputs a list of size (number of technical replicates) containing
# summary tables identical to HDX-Express output
# Behaviour of function: Performs simple automatic centroiding based on SNR
# and peak width values then offers an option for manual modification based on
# centroid inspection. The function operates interactively

# Large number of dependent functions and libraries (See to later)








mainCentMerge <- function(samplematrix,charge=1,SNR=10,bpint=0.5,sample=0,time=0){

  Samples <- unique(samplematrix[,1]) #List of Samples (E.g A1,A2...)
  NoSamples <- length(Samples)  #Number of technical replicates
  w <- c() # vector to hold widths
  cp <- c() #Vector to hold positional centroid
  c <- c() # Vector to h0ld reported centroid
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
      v <- widthFinder(p,s,bpint)
      w[i] <- (v[2]-v[1])*charge
      cp[i] <- centroidCalc(v,s)

      c[i] <- cp[i]

      interp <- linearInterp(p)

      #png(filename = paste(masses[1],"CentroidPlot","Sample",j,"-Time",TimeP[i],".png",sep=''))
      plot(s,main=paste("Sample",j,"Time",TimeP[i],sep=" "))
      plotLin(interp,p)
      plotWidth(v)
      centroidPlot(cp[i],s)
      legend('topright',c('Spectrum','Distribution Width','Centroid','Peak envelope'),
             lty=1,lwd=c(1,4,4,1),col=c('black','orange','green','red'))
      cat ("Press [Y] if satisfied with centroid, [N] if require manual editing: ")
      line <- readline()
      if(line=="Y"){

      }

      else{
        manualCentroids2(s,i,j,TimeP,masses)
      }
      #dev.off()

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





# mainCentNewMod <- function(samplematrix,charge=1,SNR=10,bpint=0.5,sample=0,time=0)
# Author: Sam MacIntyre
# Input: The feature matrix for a single peptide in the form of a large matrix
# containing all of the technical replicates
# Output: Outputs a list of size (number of technical replicates) containing
# summary tables identical to HDX-Express output
# Behaviour of function: Performs simple automatic centroiding based on SNR
# and peak width values. Saves all centroid plots into the current directory.
# Does not allow manual adjustment of peak envelopes(and consequently centroids)

# Large number of dependent functions and libraries (See to later)



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


  return(out)

}


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










importNew <- function (path = getwd(), ids = peptide.identifications,tol=0.006){
  # function to import Bruker raw spectra into r for subsequent use in HD-Express and MEMHDX
  # this function requires all spectra being copied into the current working directory (ie. getwd) in a data folder called HDXddmmyy (ie. HDX301016)
  # all spectra folders should be named according to the follwing convention:
  # sample(ie.capital letter A, B,...)_timepoint(ie. number in seconds)_replicate(ie. sequential number 1,2,3...)
  # For example:  A_0_1, A_0_2, A_0_3, .... A_3600_1, A_3600_2, A_3600_3, .... B_0_1, B_0_2, B_0_3,.... B_3600_1, B_3600_2, B_3600_3
  # Before you use this function, make sure you copy the full filepath to the clipboard (or manually enter the full filepath)
  # ie. navigate to the folder containing the raw spectra using windows explorer, right click on the file path at the top and select "Copy address as text"

  # do some prepping
  library(MALDIquant)
  library (MALDIquantForeign)
  library(reshape2)

  # define required functions

  myMALDI <- function(spectra = HDXdata_raw){
    if(all(sapply(spectra, isRegular)) != T){
      stop("Quality control issue")
    }
    spectra <- smoothIntensity(spectra,method="SavitzkyGolay",halfWindowSize=5)
    # spectra <- removeBaseline(spectra,method="SNIP", iterations=100)
    spectra <- removeBaseline(spectra,method="SNIP",iterations=100)
    #Intensity calibration (Do we need this?)
    spectra <- calibrateIntensity(spectra,method="PQN")
    # standard <- c(757.3992, 1046.6848, 1296.6848, 1347.7354, 1619.8223, 2093.0862, 2465.1983,3147.4710)
    return(spectra)
  }



  intensitymat <- function(spectra = HDXdata_processed){

    # This code extracts sample, timepoint and replicate information from spectrum metadata
    sampletable<-colsplit(string=substring(factor(sapply(spectra,function(A)metaData(A)$file)),nchar(getwd())+12,nchar(getwd())+19),"_",names=c("sample","timepoint","replicate"))
    sampletable[,3]<-substring(sampletable[,3], 1, 1)

    # This code creates the m/z vs intensity feature matrix
    mat <- matrix(,nrow = length(spectra),ncol=length(mass(spectra[[1]])))
    colnames(mat) <- mass(spectra[[1]])
    for(i in 1:length(spectra)){
      mat[i,] <- intensity(spectra[[i]])
    }
    # This script combines the two
    mat <- cbind(paste(sampletable[,1],sampletable[,3],sep=""),sampletable[,2],mat)

    return(mat)
  }

  masswindow <- function(peptide.index = 1,mat = HDXdata_matrix, masswindow = peptide.identifications[,1:2]){

    # This code extracts the m/z vs intensity values for a specific peptide masswindow from the MaldiQuant feature matrix
    min <- masswindow[peptide.index,1]
    max <- masswindow[peptide.index,2]
    flag <- 0
    for (i in 3:dim(mat)[[2]]){
      if(as.numeric(colnames(mat)[i]) > min && flag ==0){
        start <- i-1
        flag <- 1
      }
      if(as.numeric(colnames(mat)[i]) > max && flag == 1){
        end <- i
        flag <- 2
      }
    }

    # this code combines the sampleID and timepoint information with masswindow data in a new featurematrix
    # filter by sampleID and timepoint and return the ordered featurematrix
    featurematrix <- cbind(HDXdata_matrix[,1:2],(mat[,start:end]))
    featurematrix <- featurematrix[ order(featurematrix[,1], as.numeric(featurematrix[,2])), ]
    return(featurematrix)
  }

  PepWindows <- function(mat = HDXdata_matrix, mzwindow = peptide.identifications[,1:2], pepids = ids){

    # This code creates an R list object with feature matrices for each peptide mass window
    index <- c(1:length(mzwindow[,1]))
    pepwindow <- lapply(index, masswindow, mat, mzwindow)
    names(pepwindow) <- ids[,3]
    return(pepwindow)
  }

  # import all spectra
  HDXdata_raw <- importBrukerFlex(path)
  HDXdata_processed <- myMALDI()
  rm(HDXdata_raw)
  HDXdata_matrix <- intensitymat()
  rm(HDXdata_processed)
  PepID_features <- PepWindows()
  gc()
  return(PepID_features)
}








import.identifications = function (path = getwd(),filename = "Scaffold.txt"){
# function to import scaffold peptide identifications into r for subsequent use in HD-Express and MEMHDX
# Perform a Mascot database search, open the results in Scaffold, select the protein, copy and pase all peptide IDs into Notepad and save as 'Scaffold.txt'
# make sure you copy the filepath in windows explorer (or manually enter the filepath in full)
# ie. open folder containing the Scaffold.txt file in windows explorer, right click menu on top and select "Copy address as text"
# then execute the following command: PepIDs <- ScaffoldImport()

# do some prepping
library("reshape2")
library("stringr")
library ("OrgMassSpecR")
filepath = paste(path,"Scaffold.txt",sep = "/")

# define required functions

windowCalc = function(x = PepIDs)
{
# This code returns a table with the mass window (min/ max) for all identified peptides
masses <- matrix(,nrow = length(x[,1]),ncol = 2)
  for(i in 1:length(x[,1]))
  {
    masses[i,1] <- MonoisotopicMass(ConvertPeptide(x[i,1], IAA=FALSE))
    masses[i,2] <- sum(IsotopicDistributionHDX(x[i,1],incorp = 0,charge = 1, custom = list(code = "C", elements = c(C=3, H=4, N=1, O=1, S=1)))[,1][1]*(IsotopicDistributionHDX(x[i,1],incorp = 0,charge = 1,custom = list(code = "C", elements = c(C=3, H=4, N=1, O=1, S=1)))[,2]/sum(IsotopicDistributionHDX(x[i,1],incorp = 0,charge = 1,custom = list(code = "C", elements = c(C=3, H=4, N=1, O=1, S=1)))[,2])))
  }
  colnames(masses) <- c("Mono_calc","Avg_calc")
  masses<-data.frame(masses)
  return(masses)
  }

# import the Scaffold Table
PepIDs<-read.delim(filepath)


#transform some of the data
PepIDs$Sequence = substring(gsub("\\(|\\)","",PepIDs$Sequence), 2)
PepIDs$Prob = as.numeric(gsub("\\%","",PepIDs$Prob))

#remove unused columns
keeps = c("Sequence","Charge", "Observed", "Actual.Mass","Delta.PPM","Start", "Stop","NTT","Prob","P.score","Modifications","Spectrum.ID")
PepIDs = PepIDs[keeps]

# calculate the theoretical peptide masses
masses_th <- windowCalc()

# calculate the number of slow exchanging (ie.SX) and fast exchanging (ie. FX) amides
SX = nchar(PepIDs[,1])-1-as.vector(sapply(PepIDs[,1],str_count, pattern="P"))
FX = 2+as.vector(sapply(PepIDs[,1],str_count, pattern="R"))+as.vector(sapply(PepIDs[,1],str_count, pattern="K"))+as.vector(sapply(PepIDs[,1],str_count, pattern="Y"))+as.vector(sapply(PepIDs[,1],str_count, pattern="H"))+as.vector(sapply(PepIDs[,1],str_count, pattern="C"))+as.vector(sapply(PepIDs[,1],str_count, pattern="D"))+as.vector(sapply(PepIDs[,1],str_count, pattern="E"))
deltaM = PepIDs[,3] - masses_th$Mono_calc
Mono_corr = masses_th$Mono_calc + deltaM
Avg_corr = masses_th$Avg_calc + deltaM
Avg_SX = Avg_corr + (SX*1.00616)
Avg_FX = Avg_corr + ((SX+FX)*1.00616)



PepIDs <- cbind(PepIDs[,3]-3,Avg_FX+5,PepIDs,SX,FX,masses_th$Mono_calc,Mono_corr,Avg_corr,Avg_SX,Avg_FX)
colnames(PepIDs)<-c("Min","Max",keeps,"SX", "FX","Mono_calc","Mono_corr","NX_centroid","SX_centroid","FX_centroid")

gc()
return(PepIDs)
}



# Recursive function for appending lists

RecursAppend <- function(list,n){

  #basecase

  if(n==1){
    newlist <- append(list,list)
    return(newlist)
  }

  else{
    return(append(list,RecursAppend(list,n-1)))
  }
}







# takes the 3 samples for a peptide and finds
# the average value of the centroid across all timepoints

AvCent <- function(AllpeptideSamples){

  numSam <- length(AllpeptideSamples)
  numTime <- dim(AllpeptideSamples[[1]])[1]

  AvVec <- rowMeans(cbind(AllpeptideSamples[[1]][,2],AllpeptideSamples[[2]][,2],AllpeptideSamples[[3]][,2]))
  AvVec <- as.data.frame(AvVec)


  colnames(AvVec)<-"Average Centroid"

  return(AvVec)


}








