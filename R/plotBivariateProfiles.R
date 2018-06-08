# plotBivariateProfiles function
#' @author Dario Azzimonti
#' @title Plot bivariate profiles
#' @name plotBivariateProfiles
#' @description Plot bivariate profiles, for dimension up to 6.
#' @param bivProf list returned by \code{obliqueProfiles}.
#' @param allPhi a list containing the matrices Phi (dim \eqn{2xd}) for which to compute the profile extrema
#' @param Design a matrix of dimension \eqn{(2d)x numPhi} encoding the first (\code{Design[1:d,]}) and the second ((\code{Design[(d+1):(2*d),]})) axis values.
#' @param threshold if not \code{NULL} plots the level as a contour.
#' @param whichIQR which quantiles to use for the inter-quantile range plot. If \code{NULL}, automatically selects the first and the last element of \code{bivProfres_UQ$prof_quantiles_approx}
#' @param plot_options list as returned by \code{setPlotOptions}.
#' @param ... additional parameters to be passed to the plot function
#' @return plots the 2d maps of the profile sup and inf of the function for each Phi in allPhi. If threshold is not NULL also contours the threshold level.
#' @export
plotBivariateProfiles<-function(bivProf,allPhi,Design=NULL,threshold=NULL,whichIQR=NULL,plot_options=NULL,...){

  if(is.null(plot_options))
    plot_options <- setPlotOptions(plot_options = plot_options,d = bivProf$kmModel@d,num_T = 1,kmModel = bivProf$kmModel)

  trueEvals <-plot_options$fun_evals

  num_Phi<-length(allPhi)

  if(is.null(Design)){
    warning("Design not selected, auto choice not implemented")
    return(NULL)
  }

  dd_eta<-length(Design[,1])/2

  if(!is.null(threshold)){
    whichAbove = bivProf$kmModel@y>threshold
  }else{
    trueEvals = min(1,trueEvals)
  }


  ## Mean maps
  allRes <- bivProf$profMean_full

  if(plot_options$save)
    pdf(file = paste(plot_options$folderPlots,"profMean_full",plot_options$id_save,".pdf",sep=""),width = 18,height = 9)

  for(i in seq(num_Phi)){

    whichVars<-as.logical(apply(allPhi[[i]]==1,2,sum))
    pp <- bivProf$kmModel@X[,whichVars]

    if(trueEvals<3){
      zzBackgroud_max<-matrix(allRes$res$max[[i]],ncol=dd_eta)
      zzBackgroud_min<-matrix(allRes$res$min[[i]],ncol=dd_eta)
    }else{
      propPoints<-getPointProportion(pp=pp,xBins = Design[1:(dd_eta),i],
                                     yBins = Design[(dd_eta+1):(2*dd_eta),i],whichAbove = whichAbove,plt = FALSE)

      zzBackgroud_max<-propPoints$freq
      zzBackgroud_min<-propPoints$freq
    }
    layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
    image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=gray.colors(20),
          main=sprintf("Psup mean, \n Psi number %i",i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2],...)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
    switch (as.character(trueEvals),
            "0" = {},
            "1" = {points(pp,pch=3)},
            {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})


    image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=gray.colors(20),
          main=sprintf("Pinf mean, \n Psi number %i",i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2],...)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
    switch (as.character(trueEvals),
            "0" = {},
            "1" = {points(pp,pch=3)},
            {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})
  }

  if(plot_options$save)
    dev.off()


  ## quantiles maps
  if(!is.null(bivProf$res_UQ$prof_quantiles_approx)){
    for(j in seq(length(bivProf$res_UQ$prof_quantiles_approx))){
      allRes <- bivProf$res_UQ$prof_quantiles_approx[[j]]


      if(plot_options$save)
        pdf(file = paste(plot_options$folderPlots,"profQuant_",names(bivProf$res_UQ$prof_quantiles_approx)[j],plot_options$id_save,".pdf",sep=""),width = 18,height = 9)

      for(i in seq(num_Phi)){

        whichVars<-as.logical(apply(allPhi[[i]]==1,2,sum))
        pp <- bivProf$kmModel@X[,whichVars]

        if(trueEvals<3){
          zzBackgroud_max<-matrix(allRes$res$max[[i]],ncol=dd_eta)
          zzBackgroud_min<-matrix(allRes$res$min[[i]],ncol=dd_eta)
        }else{
          propPoints<-getPointProportion(pp=pp,xBins = Design[1:(dd_eta),i],
                                         yBins = Design[(dd_eta+1):(2*dd_eta),i],whichAbove = whichAbove,plt = FALSE)

          zzBackgroud_max<-propPoints$freq
          zzBackgroud_min<-propPoints$freq
        }
        layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
        image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=gray.colors(20),
              main=sprintf("Psup %s, \n Psi number %i",names(bivProf$res_UQ$prof_quantiles_approx)[j],i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2],...)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
        switch (as.character(trueEvals),
                "0" = {},
                "1" = {points(pp,pch=3)},
                {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})


        image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=gray.colors(20),
              main=sprintf("Pinf %s, \n Psi number %i",names(bivProf$res_UQ$prof_quantiles_approx)[j],i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2],...)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
        switch (as.character(trueEvals),
                "0" = {},
                "1" = {points(pp,pch=3)},
                {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})
      }

      if(plot_options$save)
        dev.off()
    }

    if(is.null(whichIQR)){
      whichIQR <- c(1,length(bivProf$res_UQ$prof_quantiles_approx))
    }

    actualQQ <- as.numeric(names(bivProf$res_UQ$prof_quantiles_approx))[whichIQR]

    # Compute weighted interquantile range
    if(!is.null(bivProf$res_UQ$prof_quantiles_approx$`0.50`$res)){
      allRes<-bivProf$res_UQ$prof_quantiles_approx$`0.50`$res
    }else{
      allRes <- bivProf$profMean_full
    }
    d<-bivProf$kmModel@d
    #/(abs(z-tt)^d+1)^(1/d)

    qqRange<-bivProf$res_UQ$prof_quantiles_approx[[whichIQR[2]]]$res
    qqRange$min<-mapply(FUN=function(x,y,z,tt){(x - y)},x=qqRange$min,y=bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$min,z=allRes$res$min,MoreArgs = list(tt=threshold))
    qqRange$max<-mapply(FUN=function(x,y,z,tt){(x - y)},x=qqRange$max,y=bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$max,z=allRes$res$max,MoreArgs = list(tt=threshold))
    #qqRange$min<-qqRange$min-bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$min
    #qqRange$max<-qqRange$max-bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$max
    ddRange<-qqRange
    ddRange$max<-ddRange$max*sapply(bivProf$res_UQ$prof_quantiles_approx[[whichIQR[2]]]$res$max,FUN=function(x){x>threshold},simplify=TRUE)*sapply(bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$max,FUN=function(x){x<threshold},simplify=TRUE)#*matrix(as.integer(qquantiles[[whichIQR[1]]]$res$max<threshold),ncol=2)

    ddRange$min<-ddRange$min*sapply(bivProf$res_UQ$prof_quantiles_approx[[whichIQR[2]]]$res$min,FUN=function(x){x>threshold},simplify=TRUE)*sapply(bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$min,FUN=function(x){x<threshold},simplify=TRUE)#matrix(as.integer(qquantiles[[whichIQR[2]]]$res$min>threshold)*as.integer(qquantiles[[whichIQR[1]]]$res$min<threshold),ncol=2)

#    ddRange$min <- ddRange$min + mapply(FUN=function(x,y,z,tt){abs(2*x - y-z)/(abs(z-tt)^d+1)^(1/d)},x=allRes$res$min,y=bivProf$res_UQ$prof_quantiles_approx[[whichIQR[2]]]$res$min,z=bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$min,MoreArgs = list(tt=threshold))
#    ddRange$max <- ddRange$max + mapply(FUN=function(x,y,z,tt){abs(2*x - y-z)/(abs(z-tt)^d+1)^(1/d)},x=allRes$res$max,y=bivProf$res_UQ$prof_quantiles_approx[[whichIQR[2]]]$res$max,z=bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$max,MoreArgs = list(tt=threshold))

    if(plot_options$save)
      pdf(file = paste(plot_options$folderPlots,"profUQrange",plot_options$id_save,".pdf",sep=""),width = 18,height = 9,pointsize = 15)

    for(i in seq(num_Phi)){

      whichVars<-as.logical(apply(allPhi[[i]]==1,2,sum))
      pp <- bivProf$kmModel@X[,whichVars]

      zzBackgroud_max<-matrix(ddRange$max[,i],ncol=dd_eta)
      zzBackgroud_min<-matrix(ddRange$min[,i],ncol=dd_eta)

      layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=gray.colors(20),
            main=sprintf("Psup mean, wIQR %.3f-%.3f, \n Psi number %i",actualQQ[1],actualQQ[2],i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2],...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})


      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=gray.colors(20),
            main=sprintf("Pinf mean, wIQR %.3f-%.3f, \n Psi number %i",actualQQ[1],actualQQ[2],i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2],...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})
    }

    if(plot_options$save)
      dev.off()



  }

  ## bound maps
  if(!is.null(bivProf$res_UQ$bound$bound)){
    for(j in seq(length(bivProf$res_UQ$bound$bound))){
      allRes <- bivProf$res_UQ$bound$bound[[j]]

      if(plot_options$save)
        pdf(file = paste(plot_options$folderPlots,"profBound_",names(bivProf$res_UQ$bound$bound)[j],plot_options$id_save,".pdf",sep=""),width = 18,height = 9)


      for(i in seq(num_Phi)){

        whichVars<-as.logical(apply(allPhi[[i]]==1,2,sum))
        pp <- bivProf$kmModel@X[,whichVars]

        if(trueEvals<3){
          zzBackgroud_max<-matrix(allRes$res$max[[i]],ncol=dd_eta)
          zzBackgroud_min<-matrix(allRes$res$min[[i]],ncol=dd_eta)
        }else{
          propPoints<-getPointProportion(pp=pp,xBins = Design[1:(dd_eta),i],
                                         yBins = Design[(dd_eta+1):(2*dd_eta),i],whichAbove = whichAbove,plt = FALSE)

          zzBackgroud_max<-propPoints$freq
          zzBackgroud_min<-propPoints$freq
        }
        layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
        image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=gray.colors(20),
              main=sprintf("Psup %s bound, \n Psi number %i",names(bivProf$res_UQ$bound$bound)[j],i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2],...)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
        switch (as.character(trueEvals),
                "0" = {},
                "1" = {points(pp,pch=3)},
                {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})


        image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=gray.colors(20),
              main=sprintf("Pinf %s bound, \n Psi number %i",names(bivProf$res_UQ$bound$bound)[j],i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2],...)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
        switch (as.character(trueEvals),
                "0" = {},
                "1" = {points(pp,pch=3)},
                {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})
      }

      if(plot_options$save)
        dev.off()
    }


    # Compute weighted bound range
    if(!is.null(bivProf$res_UQ$prof_quantiles_approx$`0.50`$res)){
      allRes<-bivProf$res_UQ$prof_quantiles_approx$`0.50`$res
    }else{
      allRes <- bivProf$profMean_full
    }
    d<-bivProf$kmModel@d
    #/(abs(z-tt)^d+1)^(1/d)

    qqRange<-bivProf$res_UQ$bound$bound$upper$res
    qqRange$min<-mapply(FUN=function(x,y,z,tt){(x - y)},x=qqRange$min,y=bivProf$res_UQ$bound$bound$lower$res$min,z=allRes$res$min,MoreArgs = list(tt=threshold))
    qqRange$max<-mapply(FUN=function(x,y,z,tt){(x - y)},x=qqRange$max,y=bivProf$res_UQ$bound$bound$lower$res$max,z=allRes$res$max,MoreArgs = list(tt=threshold))
    #qqRange$min<-qqRange$min-bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$min
    #qqRange$max<-qqRange$max-bivProf$res_UQ$prof_quantiles_approx[[whichIQR[1]]]$res$max
    ddRange<-qqRange
    ddRange$max<-ddRange$max*sapply(bivProf$res_UQ$bound$bound$upper$res$max,FUN=function(x){x>threshold},simplify=TRUE)*sapply(bivProf$res_UQ$bound$bound$lower$res$max,FUN=function(x){x<threshold},simplify=TRUE)

    ddRange$min<-ddRange$min*sapply(bivProf$res_UQ$bound$bound$upper$res$min,FUN=function(x){x>threshold},simplify=TRUE)*sapply(bivProf$res_UQ$bound$bound$lower$res$min,FUN=function(x){x<threshold},simplify=TRUE)

    #    ddRange$min <- ddRange$min + mapply(FUN=function(x,y,z,tt){abs(2*x - y-z)/(abs(z-tt)^d+1)^(1/d)},x=allRes$res$min,y=bivProf$res_UQ$bound$bound$upper$res$min,z=bivProf$res_UQ$bound$bound$lower$res$min,MoreArgs = list(tt=threshold))
    #    ddRange$max <- ddRange$max + mapply(FUN=function(x,y,z,tt){abs(2*x - y-z)/(abs(z-tt)^d+1)^(1/d)},x=allRes$res$max,y=bivProf$res_UQ$bound$bound$upper$res$max,z=bivProf$res_UQ$bound$bound$lower$res$max,MoreArgs = list(tt=threshold))

    if(plot_options$save)
      pdf(file = paste(plot_options$folderPlots,"profUBoundRange",plot_options$id_save,".pdf",sep=""),width = 18,height = 9,pointsize = 15)

    for(i in seq(num_Phi)){

      whichVars<-as.logical(apply(allPhi[[i]]==1,2,sum))
      pp <- bivProf$kmModel@X[,whichVars]

      zzBackgroud_max<-matrix(ddRange$max[,i],ncol=dd_eta)
      zzBackgroud_min<-matrix(ddRange$min[,i],ncol=dd_eta)

      layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=gray.colors(20),
            main=sprintf("Psup mean, wBound range, \n Psi number %i",i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2])#,...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})


      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=gray.colors(20),
            main=sprintf("Pinf mean, wBound range, \n Psi number %i",i),xlab= colnames(bivProf$kmModel@X)[whichVars][1], ylab = colnames(bivProf$kmModel@X)[whichVars][2])#,...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp,pch=3);points(pp[whichAbove,],pch=3,col=3)})
    }

    if(plot_options$save)
      dev.off()
  }


}
