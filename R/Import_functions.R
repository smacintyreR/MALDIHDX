


#' Import Spectra function
#' @description  Imports mass spectra and creates a list of feature matrices 
#' using the MALDIquant package
#'
#' @param path 
#' @param ids 
#' @param tol 
#' @keywords import, spectra, feature, MALDIquant
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' peptide.features <- importNew()

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



#' Imports Scaffold file for peptide identification
#' @description  Imports the Scaffold.txt file to create the an identifications
#' file 
#' @param path A singular mass spectrum
#' @param ids 
#' @param tol 
#' @keywords import, peptide, identification, Scaffold
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' PepIDs <- ScaffoldImport()

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
