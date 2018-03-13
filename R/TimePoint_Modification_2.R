# Time-point modification functions

Time_Dec_2 <- function(TimeStampVectorCSV){

  TS <- TimeStampVectorCSV[[1]]
  SplitTS <- strsplit(TS," ")
  Output <- c()
  for (i in 1:length(TS)){


    int <- SplitTS[[i]][2]

    Output[i] <- int

  }



  return(Output)

}



Times_TP_2 <- function(TimesDecimal,Init,RepNo = 3){


  minVec <- as.difftime(TimesDecimal,units="mins")
  T0 <- minVec[1]-Init

  TP <- round(minVec - T0)

  return(TP)



}



file_addTP_2 <- function(Timepoints, Sample){

  #Assumes you are in the correct working directory

  numFiles <- length(list.files())
  TP <- as.numeric(Timepoints)


  if(str_length(list.files()[1])== 1){
    for(i in 1:numFiles){

      if (i%%3 != 0){
        file.rename(as.character(i),paste(Sample,TP[i],i%%3,sep="_"))
      }

      else{
        file.rename(as.character(i),paste(Sample,TP[i],3,sep="_"))
      }


    }
  }


  else{

    for(i in 1:numFiles){
      cur <- list.files()[1]
      len <- str_length(cur)
      Split <- substring(cur,2,len)
      #split string into sample, timepoint, replicate
      file.rename(cur,paste(Sample,Split,sep=''))

    }

    i=0


  }
}




# Full script for renaming file names using appropriate conventions
# and timepoints

# Begin in Directory starting directory

RenameTP_2 <- function(ExcelSheet,initial = 4,type,sample="A") {

  timestamps <- subset(ExcelSheet,BaseName =="1r")
  timestamps <- timestamps[,"CreationTime"]



  curDate <-gsub("-","",substring(Sys.Date(),3))

  if(type == "H"){
    newFile <- paste("HDX",curDate,sep="")
  }

  else{
    newFile <- paste("DHX",curDate,sep="")
  }


  # Rename the overall folder to agree with import functions
  file.rename(list.files(),newFile)

  setwd(list.files())

  timesExtract <- Time_Dec(timestamps)

  StimesExtract <- sort(timesExtract)

  Times <- Times_TP(StimesExtract,Init = initial)

  file_addTP(Times,sample)

}




# Full script for renaming file names using appropriate conventions
# and timepoints

# Begin in Directory starting directory

RenameTPInput_2 <- function(Timepoints) {





  curDate <-gsub("-","",substring(Sys.Date(),3))


    newFile <- paste("HDX",curDate,sep="")



  # Rename the overall folder to agree with import functions
  if(identical(list.files()[1],newFile) == F ){
    file.rename(list.files()[1],newFile)
  }

  setwd(list.files()[1])

  timesExtract <- Time_Dec(timestamps)

  StimesExtract <- sort(timesExtract)

  Times <- Times_TP(StimesExtract,Init = initial)

  file_addTP(Times,Sample = sample)

  setwd("..") #Returns back to parent directory

}
