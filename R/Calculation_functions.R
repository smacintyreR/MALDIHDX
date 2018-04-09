
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



#' Peak width calculation
#' @description  Uses the spectrum peaks and a specified base peak intensity
#' to find the corresponding m/z end points
#'
#' @param spectrum A singular mass spectrum
#' @param peaks List of mass peaks for the spectrum
#' @param bpint base peak intensity (to determine width)
#' @keywords width, peaks, interpolation
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' [a,b,envelope] <- widthFinder(peaks,spectrum,bpint = 0.5)

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


#' Centroid calculation
#' @description  Calculates the centroid for a spectrum given the peak width
#'
#'
#'
#'
#' @param v Vector of end points of peak width
#' @param spectrum Mass spectrum
#' @keywords centroids,width
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' c <- centroidCalc(v,spectrum)


centroidCalc <- function(v,spectrum){

  a <- v[1]
  b <- v[2]

  new <- subset(spectrum,mass(spectrum) < b & mass(spectrum) > a)

  cent <- sum(intensity(new)*mass(new))/sum(intensity(new))

  return(cent)

}


#' Linear Equation Solver
#' @description  Finds the x coordinate of the intersection point of a linear
#' function with a specified value
#' @param p1 Vector containing first point [x1,y1]
#' @param p2 Vector containing second point [x2,y2]
#' @param value Scalar denoting required y value for intersection
#' @keywords linear, interpolation, solver, equation
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' x <- linearsolve(p1,p2,value)

linearsolve <- function(p1,p2,value){
  x1 <- p1[1]
  y1 <- p1[2]
  x2 <- p2[1]
  y2 <- p2[2]
  
  m <- (y2-y1)/(x2-x1)
  
  
  if(m != 0){


  x <- (value - y1 +m*x1)/m
  
  return(x)
  }
  
  else{
      stop("Error: gradient cannot be equal to 0")
  }

}



#' Linear Interpolation Function
#' @description  Given a set of mass peaks, finds the linear interpolation
#' between them
#' @param peaks list of mass peaks
#' @keywords linear,interpolation,peaks,
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' interp <- linearInterp(peaks)

linearInterp <- function(peaks){

  peaksmod <-cbind(mass(peaks),intensity(peaks))

  interp <- approx(x=peaksmod[,1],y=peaksmod[,2],method="linear",n=100,rule=1)
  return(interp)

}

#' Recursive Appending Function
#' @description  Given a list, recursively replicates the list for a given
#' n 
#' @param list any list
#' @param n number of times to replicate
#' @keywords replicate, append, recursive 
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' newlist <- RecursAppend(list,4)

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
