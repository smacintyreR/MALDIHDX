# Assumes that we are in the appropriate working directory (A or B, Bound or Unbound)

renameTimePoints <- function(TimepointsPositions,rep = 3){

  # Extract current folder
  dir <- getwd()
  dirlen <- str_count(dir)
  state <- substr(dir,dirlen,dirlen)

  noTimePoints <- dim(TimepointsPositions)[1]

  # Rename folders with correct formatting (E.g State_Timepoint_Replicate : A_10_2)
  for(i in 1:noTimePoints){
    setwd(list.files()[i])
    for(j in 1:rep){
      TP <- as.numeric(TimepointsPositions[i,2])
      file.rename(as.character(j),paste(state,TP,as.character(j),sep="_"))
    }
    setwd("..")
  }
}


collectFolders <- function(noReps = 3){

  # Assumes we are in the main folder containing both state folders (A and B)

  # create new directory with appropriate naming convention to store all spectra
  # i.e HDX020995

  curDate <-gsub("-","",substring(Sys.Date(),3))

  folderName <- paste("HDX",curDate,sep="")

  # Find current states (allows flexibility for more than 2 states)
  states <- list.files()
  NoStates <- length(states)

  if(dir.exists(paste("HDX",curDate,sep="")) == F){
  dir.create(folderName)
  }

  # Obtain the path of the new folder
  folderPath <- paste(getwd(),folderName,sep="/")

  # Enter A/B folders and copy correctly formatted (A_24_3) files into new folder

  for(i in 1:NoStates){
    setwd(states[i])
    subFolders <- list.files()
    noTimePoints <- length(subFolders)
    for(j in 1:noTimePoints){
      setwd(subFolders[j])
      for(k in 1:noReps){
        formattedFolders <- list.files()
        curFold <- paste(getwd(),formattedFolders[k],sep="/")
        file.copy(curFold,folderPath,recursive=T,overwrite = T)
      }
      setwd("..")
    }
    setwd("..")
  }



}


FolderDeletion <- function(){

  # Assumes we are in the main folder containing both state folders (A and B)
  origFiles <- list.files(pattern ="[A-B]")
  unlink(origFiles,recursive=T)

}
