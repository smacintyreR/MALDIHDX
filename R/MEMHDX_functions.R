# function which produces the MEMHDX formatted table for all peptides and
# and combines into a single CSV document

#GitHubTest

MEMHDXall <- function(AllpeptideCentroidsU = peptide.featuresU,
                      AllpeptideCentroidsB = peptide.featuresB,
                      Idents){

  # Extract required information from the identifications table
  IDs <- Idmodify(Idents)

  # Number of peptides to analyse
  noFeatures <- length(AllpeptideCentroidsU)

  # Calculates and stores centroids for untreated and treated sample
  CentroidsU <- lapply(AllpeptideCentroidsU,FUN = function(x) mainCentNewMod(x,SNR=3) )
  CentroidsB <- lapply(AllpeptideCentroidsB,FUN = function(x) mainCentNewMod(x,SNR=3) )





  # Builds up MEMHDX formatted matrix feature by feature - untreated
  New <- MEMHDXformat(CentroidsU[[1]],IDs = IDs,1,bound=F)
  for(i in 2:noFeatures){

    nextPep <- MEMHDXformat(CentroidsU[[i]],IDs = IDs,i,bound = FALSE)
    New <- rbind.data.frame(New,nextPep)

  }


  # Builds up MEMHDX formatted matrix feature by feature - treated
  New2 <- MEMHDXformat(CentroidsB[[1]],IDs = IDs,1)
  for(i in 2:noFeatures){

    nextPep <- MEMHDXformat(CentroidsB[[i]],IDs = IDs,i,bound = TRUE)
    New2<- rbind.data.frame(New2,nextPep)

  }

  # Combine treated and untreated matrices
  MEMHDXfull <- rbind.data.frame(New,New2)
  MEMHDXfullMat <- as.matrix(MEMHDXfull)

  # Write full table to CSV file with current date appended
  write.csv(MEMHDXfullMat,file=paste("MEMHDX",Sys.Date(),".csv",sep=""))
  return(MEMHDXfull)

}


# function which combines the peptide.identification information and centroids
# for one peptide

MEMHDXformat <- function(Centroids,IDs,pepno,bound = TRUE){


  # Number of technical replicates
  noRep <- length(Centroids)

  # Number of Timepoints
  noTime <- dim(Centroids[[1]])[1]


  # Create 0 matrix with appropriate dimensions
  CentroidMatrix <- matrix(0,noTime*noRep,5)
  CentroidMatrix <- as.data.frame(CentroidMatrix)

  k = 0
  for(i in 1:noTime){
    for(j  in 1:noRep){
      k = k+1
      # Fill MEMHDX style matrix
      CentroidMatrix[k,4] <- j
      CentroidMatrix[k,5] <- Centroids[[j]][i,2]
      CentroidMatrix[k,2] <- as.character(Centroids[[j]][i,1])

    }

  }


  # Make this more flexible to accommodate more states?
  if(bound == T){
    state <- "BOUND"
  }

  else{
    state <- "UNBOUND"
  }

  stateVec <- rep.row(state, noTime*noRep)


  # Replicate ID information
  IDsPep <- IDs[pepno,]
  namesID <-colnames(IDsPep)
  IDsPep <- rep.row(IDsPep,noTime*noRep)

  # Extract charge column
  charge <- IDsPep[,5]
  IDsPep <- IDsPep[,-5]

  # Create final matrix
  IDsPep <- as.data.frame(IDsPep)
  colnames(IDsPep) <- namesID[-5]
  CentroidMatrix[,1] <- stateVec
  CentroidMatrix[,3] <- charge

  # Set column names
  colnames(CentroidMatrix) <- c("State","Exposure","z","Replicate","Center")



  # Combine identifications and centroids
  FullMatrix <- cbind(IDsPep,CentroidMatrix)

  return(FullMatrix)

}







# function to modify peptide.identification table
# Extract columns required to format MEMHDX table

Idmodify <- function(IDs = peptide.identifications){

  # Subsets necessary columns for MEMHDX table
  New <- IDs[,c(8,9,3,15,4)]

  # Modify colnames to adhere to MEMHDX format convention
  colnames(New)[2] <- "End"
  colnames(New)[4] <- "MaxUptake"

  return(New)

}
