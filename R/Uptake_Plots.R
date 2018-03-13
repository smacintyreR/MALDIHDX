

#' Linear regression and plotting
#'
#' @description Calculates linear regression coefficients for HDX and DHX
#' uptakes, plots the curves and returns the regression coefficients for
#' the two curves
#' @param pepnum
#' @param pepIDS
#' @param HDXcents
#' @param DXcents
#' @param times
#' @keywords linear, regression, uptakes
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' coefs <- PlotDHLinear()
PlotDHLinear <- function(pepnum,pepIDS = peptide.identifications,HDXcents,DXcents,times){

  # HDX Uptake values
  UptakesHDX <- (HDXcents - pepIDS[pepnum,19])/(pepIDS[pepnum,20]-pepIDS[pepnum,19])*100
  UptakesHDX <- as.data.frame(UptakesHDX)
  # DX uptake values
  UptakesDX <- (DXcents - pepIDS[pepnum,19])/(pepIDS[pepnum,21]-pepIDS[pepnum,19])*100
  UptakesDX <- as.data.frame(UptakesDX)

  # Plotting
  plot(times,UptakesHDX[,1],ylim=c(0,100),xlim=c(0,30),col = "blue",pch=4,ylab="Deuterium Uptake (%)",xlab = "Time (minutes)",main = "DHX vs HDX (Linear)")
  points(times,UptakesDX[,1],col="red",pch=3)
  legend('bottomright',c('HDX','DX'),
         pch=c(4,3),col=c('blue','red'))


  # Calculation of linear
  ModelD <- lm(UptakesDX[,1]~times)
  Model <- lm(UptakesHDX[,1]~times)

  # DX Trendline
  abline(ModelD,col='red',lty=2)
  # HDX Trendline
  abline(Model,col='blue',lty=2)

  # Return coefficients
  return(rbind(coef(Model),coef(ModelD)))

}



#' Polynomial (order 2) regression and plotting
#'
#' @description Calculates polynomial (order 2) regression coefficients for HDX
#'  and DHX uptakes, plots the curves and returns the regression coefficients for
#' the two curves
#' @param pepnum
#' @param pepIDS
#' @param HDXcents
#' @param DXcents
#' @param times
#' @keywords polynomial, regression, uptakes
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' coefs <- PlotDHPoly2()

PlotDHPoly2 <- function(pepnum,pepIDS = peptide.identifications,HDXcents,DXcents,times){

  # HDX Uptake values
  UptakesHDX <- (HDXcents - pepIDS[pepnum,19])/(pepIDS[pepnum,20]-pepIDS[pepnum,19])*100
  UptakesHDX <- as.data.frame(UptakesHDX)
  # DX uptake values
  UptakesDX <- (DXcents - pepIDS[pepnum,19])/(pepIDS[pepnum,21]-pepIDS[pepnum,19])*100
  UptakesDX <- as.data.frame(UptakesDX)


  # Plotting
  plot(times,UptakesHDX[,1],ylim=c(0,100),xlim=c(0,30),col = "blue",pch=4,ylab="Deuterium Uptake (%)",xlab = "Time (minutes)",main = "DHX vs HDX (Polynomial Order 2)")
  points(times,UptakesDX[,1],col="red",pch=3)
  legend('bottomright',c('HDX','DX'),
         pch=c(4,3),col=c('blue','red'))

  # General Polynomial
  f <- function(x,a,b,d) {(a*x^2) + (b*x) + d}

  # HDX Trendline
  fitH <- nls(UptakesHDX[,1] ~ f(times,a,b,d), start = c(a=1, b=1, d=1))
  coH <- coef(fitH)
  curve(f(x, a=coH[1], b=coH[2], d=coH[3]), add = TRUE, col="green", lwd=2)

  # DX Trendline
  fitD <- nls(UptakesDX[,1] ~ f(times,a,b,d), start = c(a=1, b=1, d=1))
  coD <- coef(fitD)
  curve(f(x, a=coD[1], b=coD[2], d=coD[3]), add = TRUE, col="pink", lwd=2)

  # Return coefficients
  return(rbind(coH,coD))

}



#' Logarithmic regression and plotting
#'
#' @description Calculates logarithmic regression coefficients for HDX
#'  and DHX
#' uptakes, plots the curves and returns the regression coefficients for
#' the two curves - uses equation of the form: (a*log(x) + b)
#' @param pepnum
#' @param pepIDS
#' @param HDXcents
#' @param DXcents
#' @param times
#' @keywords polynomial, regression, uptakes
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' coefs <- PlotDHLog()

PlotDHLog <- function(pepnum,pepIDS = peptide.identifications,HDXcents,DXcents,times){

  # HDX Uptake values
  UptakesHDX <- (HDXcents - pepIDS[pepnum,19])/(pepIDS[pepnum,20]-pepIDS[pepnum,19])*100
  UptakesHDX <- as.data.frame(UptakesHDX)
  # DX uptake values
  UptakesDX <- (DXcents - pepIDS[pepnum,19])/(pepIDS[pepnum,21]-pepIDS[pepnum,19])*100
  UptakesDX <- as.data.frame(UptakesDX)

  # Plotting
  plot(times,UptakesHDX[,1],ylim=c(0,100),xlim=c(0,30),col = "blue",pch=4,ylab="Deuterium Uptake (%)",xlab = "Time (minutes)",main = "DHX vs HDX (Log)")
  points(times,UptakesDX[,1],col="red",pch=3)
  legend('bottomright',c('HDX','DX'),
         pch=c(4,3),col=c('blue','red'))

  # General Log Function
  f <- function(x,a,b) {a * log(x) + b}

  # HDX Trendline
  fitH <- nls(UptakesHDX[,1] ~ f(times,a,b), start = c(a=1, b=1))
  coH <- coef(fitH)
  # Plotting
  curve(f(x, a=coH[1], b=coH[2]), add = TRUE, col="orange", lwd=2)

  # DX Trendline
  fitD <- nls(UptakesDX[,1] ~ f(times,a,b), start = c(a=1, b=1))
  coD <- coef(fitD)
  # Plotting
  curve(f(x, a=coD[1], b=coD[2]), add = TRUE, col="purple", lwd=2)

  # Return coefficients
  return(rbind(coH,coD))

}



#' Logarithmic intersection
#'
#' @description Calculates the intersection point of two logarithmic functions
#' given their coefficients
#' @param pepnum
#' @param pepIDS
#' @param HDXcents
#' @param DXcents
#' @param times
#' @keywords log,logarithmic, intersection
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' [x,y] <- LogIntersect(LogCoefs)

LogIntersect <- function(LogCoefs){

  # Calculation of intersection
  x <- exp((LogCoefs[2,2]-LogCoefs[1,2])/(LogCoefs[1,1]-LogCoefs[2,1]))
  y <- LogCoefs[1,1]*log(x) + LogCoefs[1,2]


  intersection <- c(x,y)

  # Check x is between 0.001 and 30 and is a real value
  if(x > 0.001 & x <=30 & !is.na(x)){

    return(intersection)
  }

  else{
    return(FALSE)
  }
}


#' Polynomial intersection
#'
#' @description Calculates the intersection point of two polynomial functions
#' given their coefficients
#' @param pepnum
#' @param pepIDS
#' @param HDXcents
#' @param DXcents
#' @param times
#' @keywords log,logarithmic, intersection
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' [x,y] <- PolyIntersect(PolyCoefs)

PolyIntersect <- function(PolyCoefs){

  a <- PolyCoefs[1,1] - PolyCoefs[2,1]
  b <- PolyCoefs[1,2] - PolyCoefs[2,2]
  c <- PolyCoefs[1,3] - PolyCoefs[2,3]

  # Calculation of intersection
  x <- (-b - sqrt(b^2 - 4*a*c))/(2*a)
  y <- PolyCoefs[1,1]*x^2 + PolyCoefs[1,2]*x + PolyCoefs[1,3]

  intersection <- c(x,y)

  # Check x is between 0 and 30 and is a real value
  if(x >=0 & x <=30 & !is.na(x)){

    return(intersection)
  }

  else{
    return(FALSE)
  }


}



#' Linear intersection
#'
#' @description Calculates the intersection point of two linear functions
#' given their coefficients
#' @param LinCoefs
#' @keywords linear, intersection
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' [x,y] <- LinIntersect(LinCoefs)

LinIntersect <- function(LinCoefs){

  x <- (LinCoefs[2,1] - LinCoefs[1,1])/(LinCoefs[1,2] - LinCoefs[2,2])
  y <- LinCoefs[1,1] + LinCoefs[1,2]*x

  intersection <- c(x,y)

  if(x >=0 & x <=30 & !is.na(x)){

    return(intersection)
  }

  else{
    return(FALSE)
  }

}

#' Average centroids for a complete HDX/DHX experiment
#'
#' @description Calculates average centroids across each sample and each time
#' point
#' @param LinCoefs
#' @keywords linear, intersection
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' [x,y] <- LinIntersect(LinCoefs)


AverageCents <- function(PepDX=peptide.featuresD,PepHDX = peptide.features,SNR = 5,pepnums){

  CentsHD <- c()
  CentD <- c()
  for(i in 1:pepnums){
    CentsHD[i] <- AvCent(mainCentNewMod(PepHDX[[i]],SNR=SNR,bpint=0.3))
    CentD[i] <- AvCent(mainCentNewMod(PepDX[[i]],SNR=SNR,bpint=0.3))

  }

  Allcents <- list(CentD,CentsHD)

  return(Allcents)

}




#' All regression plots for a single peptide
#'
#' @description Calculates the regression curves (linear, logarithmic and
#' polynomial) for a full peptide
#'
#' @param pepnum
#' @param pepIDs
#' @param HDXcents
#' @param DXcents
#' @param times
#' @keywords linear, regression, polynomial, logarithmic, pdf, plots
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' PlotDHAll(pepnum,DXcents,DXcents,times)



PlotDHAll <- function(pepnum,pepIDS = peptide.identifications,HDXcents,DXcents,times){
  pdf(paste("DXvsHDXRegressionPlotsPeptide",pepnum,".pdf"))
  old.par <- par(mfrow=c(2,2))

  LinCoefs <- PlotDHLinear(pepnum,pepIDS,HDXcents,DXcents,times)
  LinInt <- LinIntersect(LinCoefs)
  LogCoefs <- PlotDHLog(pepnum,pepIDS,HDXcents,DXcents,times)
  LogInt <- LogIntersect(LogCoefs)
  PolyCoefs <- PlotDHPoly2(pepnum,pepIDS,HDXcents,DXcents,times)
  PolyInt <- PolyIntersect(PolyCoefs)
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

  if(LinInt){
  text(0.5,0.8,paste("Linear Intersection: x = ",round(LinInt[1],2),"y = ",round(LinInt[2],2)),cex=0.9)
  }
  if(LogInt){
  text(0.5,0.6,paste("Log Intersection: x = ",round(LogInt[1],2),"y = ",round(LogInt[2],2)),cex=0.9)
  }
  if(PolyInt){
  text(0.5,0.4,paste("Polynomial Intersection: x = ",round(PolyInt[1],2),"y = ",round(PolyInt[2],2)),cex=0.9)
  }
  dev.off()
  par(old.par)

}


#' All regression plots for a full HDX/DHX experiment
#'
#' @description Calculates the regression curves (linear, logarithmic and
#' polynomial) for a complete HDX/DHX experiment and saves each set of 3 plots
#' in a pdf. All pdfs are saved in a directory called: "AllRegressionPlots"
#' @param Allcents
#' @param pepIDs
#' @param pepnumber
#' @param times
#' @keywords linear, regression, polynomial, logarithmic, pdf, plots
#' @author Sam MacIntyre \email{smacintyre@@csiro.au}
#' @export
#' @examples
#' AllRegressionPlots(Allcents,pepnumber,times)


AllRegressionPlots <- function(Allcents,pepIDs = peptide.identifications,pepnumber,times){



  CentsD <- Allcents[[1]]
  CentsHD <- Allcents[[2]]

  dir.create("AllRegressionPlots")
  setwd("AllRegressionPlots")
  for(i in 1:pepnumber){

    PlotDHAll(i,pepIDS = pepIDs,CentsHD[[i]],CentsD[[i]],times)


  }


}











