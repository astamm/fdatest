#' Plot the p-values of the ITP1 test
#' 
#' @param x An object of class `ITP1`.
#' @param xrange Range of the x-axis.
#' @param alpha1 Significance level for the first type of error.
#' @param alpha2 Significance level for the second type of error.
#' @param ylab Label for the y-axis.
#' @param main Title of the plot.
#' @param lwd Line width.
#' @param col Color of the lines.
#' @param pch Symbol for the points.
#' @param ylim Range of the y-axis.
#' @param ... Additional arguments to be passed to the `matplot` function.
#' 
#' @return A plot of the p-values of the ITP1 test.
#' 
#' @export
plot.ITP1 <- function(x, 
                      xrange = c(0, 1),
                      alpha1 = 0.05, 
                      alpha2 = 0.01,
                      ylab = "Functional Data", 
                      main = NULL, 
                      lwd = 1, 
                      col = 1, 
                      pch = 16, 
                      ylim = range(object$data.eval),
                      ...) {
  if(alpha1 < alpha2){
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  object <- x
  
  graphics::par(ask=TRUE) 
  if(object$basis=='Fourier'){
    p <- length(object$pval)
    J <- dim(object$data.eval)[2]
    n <- dim(object$data.eval)[1]
    xmin <- xrange[1]
    xmax <- xrange[2]
    abscissa.pval = 1:p
    Abscissa = seq(xmin,xmax,len=J)
    main.data <- paste(main,': Functional Data')
    main.data <- sub("^ : +", "", main.data)
    fda::matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=col,lwd=lwd,ylim=ylim,...)
    if(length(object$mu)==1){
      abscissa.mu <- Abscissa
      mu <- rep(object$mu,1000)
    }else{
      Abscissa <- seq(xmin,xmax,length.out=length(object$mu))
      mu <- object$mu
    }
    graphics::lines(abscissa.mu,mu,col='gray',lwd=2)
    
    ################################################################
    # pval
    main.p <- paste(main,': Adjusted p-values')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$adjusted.pval,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',xlab='Frequency')
    difference1 <- which(object$adjusted.pval<alpha1)
    if(length(difference1)>0){
      for(j in 1:length(difference1)){
        min.rect <- abscissa.pval[difference1[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    difference2 <- which(object$adjusted.pval<alpha2)
    if(length(difference2)>0){
      for(j in 1:length(difference2)){
        min.rect <- abscissa.pval[difference2[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray80',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(1:p,object$adjusted.pval,pch=pch)
    
    
    
  }else if(object$basis=='B-spline'){  
    p <- length(object$pval)
    J <- dim(object$data.eval)[2]
    n <- dim(object$data.eval)[1]
    xmin <- xrange[1]
    xmax <- xrange[2]
    abscissa.pval = seq(xmin,xmax,len=p)
    Abscissa = seq(xmin,xmax,len=J)
    main.data <- paste(main,': Functional Data')
    main.data <- sub("^ : +", "", main.data)
    fda::matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=col,lwd=lwd,ylim=ylim,...)
    difference1 <- which(object$adjusted.pval<alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$adjusted.pval<alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray80", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    fda::matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=col,lwd=lwd,add=TRUE,...)
    if(length(object$mu)==1){
      abscissa.mu <- Abscissa
      mu <- rep(object$mu,1000)
    }else{
      Abscissa <- seq(xmin,xmax,length.out=length(object$mu))
      mu <- object$mu
    }
    graphics::lines(abscissa.mu,mu,col='gray',lwd=2)
    
    ################################################################
    # pval
    main.p <- paste(main,': Adjusted p-values')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$adjusted.pval,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',...)
    difference1 <- which(object$adjusted.pval<alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$adjusted.pval<alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray80", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(abscissa.pval,object$adjusted.pval,pch=pch)
    
  }
  graphics::par(ask=FALSE) 
}
