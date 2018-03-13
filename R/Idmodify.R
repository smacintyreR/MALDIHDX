# function to modify peptide.identification table
# Extract columns required to format MEMHDX table

Idmodify <- function(IDs = peptide.identifications){

  New <- IDs[,c(8,9,3,15,4)]

  return(New)

}




















