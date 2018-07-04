


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

importNew <- function (path =getwd(), ids = peptide.identifications,tol=0.006){
    # function to import Bruker raw spectra into r for subsequent use in HD-Express and MEMHDX
    # this function requires all spectra being copied into the current working directory (ie. getwd) in a data folder called HDXddmmyy (ie. HDX301016)
    # all spectra folders should be named according to the follwing convention:
    # sample(ie.capital letter A, B,...)_timepoint(ie. number in seconds)_replicate(ie. sequential number 1,2,3...)
    # For example:  A_0_1, A_0_2, A_0_3, .... A_3600_1, A_3600_2, A_3600_3, .... B_0_1, B_0_2, B_0_3,.... B_3600_1, B_3600_2, B_3600_3
    # Before you use this function, make sure you copy the full filepath to the clipboard (or manually enter the full filepath)
    # ie. navigate to the folder containing the raw spectra using windows explorer, right click on the file path at the top and select "Copy address as text"
    
    
    # import all spectra
    HDXdata_raw <- importBrukerFlex(path)
    HDXdata_processed <- myMALDI(HDXdata_raw)
    HDXdata_matrix <- intensitymat(HDXdata_processed)
    PepID_features <- PepWindows(HDXdata_matrix,ids)
    gc()
    return(PepID_features)
}



PepWindows <- function(HDXdata_matrix,Idents = ids){
    
    mzwindow <- Idents[,1:2]
    
    # This code creates an R list object with feature matrices for each peptide mass window
    index <- c(1:length(mzwindow[,1]))
    pepwindow <- lapply(index, masswindow , HDXdata_matrix, mzwindow)
    names(pepwindow) <- Idents[,4  ]
    return(pepwindow)
}





myMALDI <- function(spectra){
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



masswindow <- function(peptide.index = 1,HDXdata_matrix, masswindow = ids[,1:2]){
    
    # This code extracts the m/z vs intensity values for a specific peptide masswindow from the MaldiQuant feature matrix
    min <- masswindow[peptide.index,1]
    max <- masswindow[peptide.index,2]
    flag <- 0
    
    for (i in 3:dim(HDXdata_matrix)[2]){
        
        curMass <- as.numeric(colnames(HDXdata_matrix)[i])
        
    
        if(curMass > min && flag ==0){
            print(curMass > min)
            start <- i-1
            flag <- 1
            
         
        }
        if(curMass > max && flag == 1){
            end <- i
            flag <- 2
        }
    }

    
    # this code combines the sampleID and timepoint information with masswindow data in a new featurematrix
    # filter by sampleID and timepoint and return the ordered featurematrix
    featurematrix <- cbind(HDXdata_matrix[,1:2],(HDXdata_matrix[,start:end]))
    featurematrix <- featurematrix[ order(featurematrix[,1], as.numeric(featurematrix[,2])), ]
    return(featurematrix)
}   




intensitymat <- function(spectra){
    
    # This code extracts sample, timepoint and replicate information from spectrum metadata
    sampletable<-colsplit(string=str_extract(factor(sapply(spectra,function(A)metaData(A)$file)),"[A-Z]_\\d+\\d*_\\d+\\d*"),"_",names=c("sample","timepoint","replicate"))
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


  filepath = paste(path,"Scaffold.csv",sep = "/")

  # define required functions

 

  # import the Scaffold Table
  PepIDs<-read.csv(filepath)

PepIDs[,1] <- as.character(PepIDs[,1])


  # calculate the number of slow exchanging (ie.SX) and fast exchanging (ie. FX) amides
  SX = nchar(PepIDs[,1])-1-as.vector(sapply(PepIDs[,1],str_count, pattern="P"))




  PepIDs <- cbind(PepIDs[,2]-3,PepIDs[,2] +SX +PepIDs[,2]*0.003,SX,PepIDs)
  colnames(PepIDs)[c(1,2,3)]<-c("Min","Max","MaxUptake")
  return(PepIDs)
}


DFtoSpec <- function(peptideTable){
    
    samples <- peptideTable[,1]
    TP <- peptideTable[,2]
    
    new <- peptideTable[,-c(1,2)]
    spectra <- list()
    
    Index <- dim(peptideTable)[1]
    
    for(i in 1:Index){
        
        spectra[[i]] <- createMassSpectrum(as.numeric(colnames(new)),as.numeric(new[i,]),metaData = list(StateRep = samples[i],time=TP[i]))
        
    }
    
    return(spectra)
    
    
    
}
