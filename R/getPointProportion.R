# getPointProportion function
#' @author Dario Azzimonti
#' @title Obtain proportion of true observations in excursion set
#' @name getPointProportion
#' @description Computes the proportion of observations in the excursion set from true function evaluations,
#'  binned by the grid determined with \code{xBins}, \code{yBins}.
#' @param pp a matrix of dimension nPts x 2 with the true points locations in 2 dimensions. 
#' @param xBins numerical vector with the ordered breaks of the grid along the x axis
#' @param yBins numerical vector with the ordered breaks of the grid along the y axis
#' @param whichAbove boolean vector of dimension nPts, selects the points above
#' @param plt if not \code{TRUE} plots the grid, the points and the counts for each cell.
#' @return a list containing \code{above}, the counts of points in excursion, \code{full} the counts per cell of all points,
#'  \code{freq}, the relative frequence. 
#' @export
getPointProportion <- function(pp,xBins, yBins, whichAbove, plt=FALSE){
  
  nnR <- length(xBins)
  nnC <- length(yBins)
  
  # get bins for full pp 
  binxy <- data.frame(x=findInterval(pp[,1], xBins),
                      y=findInterval(pp[,2], yBins))
  
  init_table<- table(binxy)
  full_table<-matrix(0,ncol=(nnC-1),nrow=(nnR-1))
  full_table[sort(unique(binxy$x)),sort(unique(binxy$y))]<-init_table
  colnames(full_table) <- 1:(nnC-1)
  rownames(full_table) <- 1:(nnR-1)
  
  if(plt){
    plot(expand.grid(xBins,yBins), t="n", xaxs="i", yaxs="i",main="full")
    points(pp, col="blue", pch="+") 
    abline(v=xBins, h=yBins)
    
    ddd <- as.data.frame.table(full_table)
    names(ddd)[1:2]<-c("x","y")
    xxx <- xBins[-length(xBins)] + 0.5*diff(xBins)
    ddd$x <- xxx[ddd$x]
    yyy <- yBins[-length(yBins)] + 0.5*diff(yBins)
    ddd$y <- yyy[ddd$y]
    with(ddd, text(x, y, label=Freq))
  }
  
  # get bins for above pp 
  binxy <- data.frame(x=findInterval(pp[whichAbove,1], xBins),
                      y=findInterval(pp[whichAbove,2], yBins))
  
  init_table<- table(binxy)
  above_table<-matrix(0,ncol=(nnC-1),nrow=(nnR-1))
  above_table[sort(unique(binxy$x)),sort(unique(binxy$y))]<-init_table
  colnames(above_table) <- 1:(nnC-1)
  rownames(above_table) <- 1:(nnR-1)
  
  freqEx<-above_table/full_table
  freqEx[is.nan(freqEx)]<-rep(0,sum(is.nan(freqEx)))
  
  if(plt){
    plot(expand.grid(xBins,yBins), t="n", xaxs="i", yaxs="i",main="whichAbove")
    points(pp[whichAbove,], col="blue", pch="+") 
    abline(v=xBins, h=yBins)
    
    ddd <- as.data.frame.table(above_table)
    names(ddd)[1:2]<-c("x","y")
    xxx <- xBins[-length(xBins)] + 0.5*diff(xBins)
    ddd$x <- xxx[ddd$x]
    yyy <- yBins[-length(yBins)] + 0.5*diff(yBins)
    ddd$y <- yyy[ddd$y]
    with(ddd, text(x, y, label=Freq))
  }
  
  return(list(above=above_table,full=full_table,freq=freqEx))
  
}