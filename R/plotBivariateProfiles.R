# plotBivariateProfiles function
#' @author Dario Azzimonti
#' @title Plot bivariate profiles
#' @name plotBivariateProfiles
#' @description Plot bivariate profiles, for dimension up to 6.
#' @param bivProf list returned by \code{obliqueProfiles}.
#' @param allPsi a list containing the matrices Psi (dim \eqn{2xd}) for which to compute the profile extrema
#' @param Design a matrix of dimension \eqn{(2d)x numPsi} encoding the first (\code{Design[1:d,]}) and the second ((\code{Design[(d+1):(2*d),]})) axis values.
#' @param threshold if not \code{NULL} plots the level as a contour.
#' @param whichIQR which quantiles to use for the inter-quantile range plot. If \code{NULL}, automatically selects the first and the last element of \code{bivProfres_UQ$prof_quantiles_approx}
#' @param plot_options list as returned by \code{setPlotOptions}.
#' @param ... additional parameters to be passed to the plot function
#' @return plots the 2d maps of the profile sup and inf of the function for each Psi in allPsi. If threshold is not NULL also contours the threshold level.
#' @export
plotBivariateProfiles<-function(bivProf,allPsi,Design=NULL,threshold=NULL,whichIQR=NULL,plot_options=NULL,...){

  if(is.null(plot_options))
    plot_options <- setPlotOptions(plot_options = plot_options,d = bivProf$kmModel@d,num_T = 1,kmModel = bivProf$kmModel)

  trueEvals <-plot_options$fun_evals

  num_Psi<-length(allPsi)

  if(is.null(Design)){
    warning("Design not selected, auto choice not implemented")
    return(NULL)
  }

  dd_eta<-length(Design[,1])/2

  bkg_col <- gray.colors(20)

  if(!is.null(threshold)){
    whichAbove = bivProf$kmModel@y>threshold
  }else{
    trueEvals = min(1,trueEvals)
  }


  ## Mean maps
  allRes <- bivProf$profMean_full

  if(plot_options$save)
    pdf(file = paste(plot_options$folderPlots,"profMean_full",plot_options$id_save,".pdf",sep=""),width = 18,height = 9)

  for(i in seq(num_Psi)){

#    whichVars<-as.logical(apply(allPsi[[i]]==1,2,sum))
#    pp <- bivProf$kmModel@X[,whichVars]

    wwVV<- apply(allPsi[[i]],2,function(x)return(x>0))

    pp<-bivProf$kmModel@X%*%t(allPsi[[i]])

    xlab_i <- paste(allPsi[[i]][1,wwVV[1,]],"*",plot_options$coord_names[wwVV[1,]],collapse = " + ")
    ylab_i <- paste(allPsi[[i]][2,wwVV[2,]],"*",plot_options$coord_names[wwVV[2,]],collapse = " + ")

    if(trueEvals<3){
      zzBackgroud_max<-matrix(allRes$res$max[[i]],ncol=dd_eta)
      zzBackgroud_min<-matrix(allRes$res$min[[i]],ncol=dd_eta)
      bck_name <- "mean"
    }else{
      propPoints<-getPointProportion(pp=pp,xBins = Design[1:(dd_eta),i],
                                     yBins = Design[(dd_eta+1):(2*dd_eta),i],whichAbove = whichAbove,plt = FALSE)

      zzBackgroud_max<-propPoints$freq
      zzBackgroud_min<-propPoints$freq
      bck_name <- "frequency"
    }
#    layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
    layout(matrix(1:4,nrow=1),widths=c(0.45,0.05,0.45,0.05))
    par(mar=c(4,4,4,0.2) + 0.2)
    image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=bkg_col,
          main=bquote(P[Psi[.(i)]]^sup ~" mean, "~  .(bck_name) ),xlab= xlab_i, ylab = ylab_i,...)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
    switch (as.character(trueEvals),
            "0" = {},
            "1" = {points(pp,pch=3)},
            {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})

    legend_info<-hist(zzBackgroud_max,plot=FALSE,breaks=20)
    which_col_leg <- which(legend_info$counts>0)
    par(mar=c(5,0.2,5,2.8))
    image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
    axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))

    par(mar=c(4,4,4,0.2) + 0.2)
    image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=bkg_col,
          main=bquote(P[Psi[.(i)]]^inf ~" mean, "~  .(bck_name) ),xlab= xlab_i, ylab = ylab_i,...)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
    switch (as.character(trueEvals),
            "0" = {},
            "1" = {points(pp,pch=3)},
            {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})

    legend_info<-hist(zzBackgroud_min,plot=FALSE,breaks=20)
    which_col_leg <- which(legend_info$counts>0)
    par(mar=c(5,0.2,5,2.8))
    image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
    axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))

  }

  if(plot_options$save)
    dev.off()


  ## quantiles maps
  if(!is.null(bivProf$res_UQ$prof_quantiles_approx)){
    for(j in seq(length(bivProf$res_UQ$prof_quantiles_approx))){
      allRes <- bivProf$res_UQ$prof_quantiles_approx[[j]]


      if(plot_options$save)
        pdf(file = paste(plot_options$folderPlots,"profQuant_",names(bivProf$res_UQ$prof_quantiles_approx)[j],plot_options$id_save,".pdf",sep=""),width = 18,height = 9)

      for(i in seq(num_Psi)){

#        whichVars<-as.logical(apply(allPsi[[i]]==1,2,sum))
#        pp <- bivProf$kmModel@X[,whichVars]

        wwVV<- apply(allPsi[[i]],2,function(x)return(x>0))

        pp<-bivProf$kmModel@X%*%t(allPsi[[i]])

        xlab_i <- paste(allPsi[[i]][1,wwVV[1,]],"*",plot_options$coord_names[wwVV[1,]],collapse = " + ")
        ylab_i <- paste(allPsi[[i]][2,wwVV[2,]],"*",plot_options$coord_names[wwVV[2,]],collapse = " + ")

        if(trueEvals<3){
          zzBackgroud_max<-matrix(allRes$res$max[[i]],ncol=dd_eta)
          zzBackgroud_min<-matrix(allRes$res$min[[i]],ncol=dd_eta)
          bck_name <- "mean"
        }else{
          propPoints<-getPointProportion(pp=pp,xBins = Design[1:(dd_eta),i],
                                         yBins = Design[(dd_eta+1):(2*dd_eta),i],whichAbove = whichAbove,plt = FALSE)

          zzBackgroud_max<-propPoints$freq
          zzBackgroud_min<-propPoints$freq
          bck_name <- "frequency"
        }
#        layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
        layout(matrix(1:4,nrow=1),widths=c(0.45,0.05,0.45,0.05))
        par(mar=c(4,4,4,0.2) + 0.2)
        image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=bkg_col,
              main=bquote(P[Psi[.(i)]]^sup ~" Z (" ~ .(names(bivProf$res_UQ$prof_quantiles_approx)[j]) ~ "), "~  .(bck_name) ),xlab= xlab_i, ylab = ylab_i,...)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
        switch (as.character(trueEvals),
                "0" = {},
                "1" = {points(pp,pch=3)},
                {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})


        legend_info<-hist(zzBackgroud_max,plot=FALSE,breaks=20)
        which_col_leg <- which(legend_info$counts>0)
        par(mar=c(5,0.2,5,2.8))
        image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
        axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))

        par(mar=c(4,4,4,0.2) + 0.2)
        image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=bkg_col,
              main=bquote(P[Psi[.(i)]]^inf ~" Z (" ~ .(names(bivProf$res_UQ$prof_quantiles_approx)[j]) ~ "), " ~ .(bck_name) ),xlab= xlab_i, ylab = ylab_i,...)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
        switch (as.character(trueEvals),
                "0" = {},
                "1" = {points(pp,pch=3)},
                {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})


        legend_info<-hist(zzBackgroud_min,plot=FALSE,breaks=20)
        which_col_leg <- which(legend_info$counts>0)
        par(mar=c(5,0.2,5,2.8))
        image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
        axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg))#,mgp=c(0,.5,0))
      }

      if(plot_options$save)
        dev.off()
    }

    if(is.null(whichIQR)){
      whichIQR <- c(1,length(bivProf$res_UQ$prof_quantiles_approx))
    }

    actualQQ <- round(as.numeric(names(bivProf$res_UQ$prof_quantiles_approx))[whichIQR],3)

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

    for(i in seq(num_Psi)){

#      whichVars<-as.logical(apply(allPsi[[i]]==1,2,sum))
#      pp <- bivProf$kmModel@X[,whichVars]

      wwVV<- apply(allPsi[[i]],2,function(x)return(x>0))

      pp<-bivProf$kmModel@X%*%t(allPsi[[i]])

      xlab_i <- paste(allPsi[[i]][1,wwVV[1,]],"*",plot_options$coord_names[wwVV[1,]],collapse = " + ")
      ylab_i <- paste(allPsi[[i]][2,wwVV[2,]],"*",plot_options$coord_names[wwVV[2,]],collapse = " + ")

      zzBackgroud_max<-matrix(ddRange$max[,i],ncol=dd_eta)
      zzBackgroud_min<-matrix(ddRange$min[,i],ncol=dd_eta)

#      layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
      layout(matrix(1:4,nrow=1),widths=c(0.45,0.05,0.45,0.05))
      par(mar=c(4,4,4,0.2) + 0.2)
      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=bkg_col,
            main=bquote(P[Psi[.(i)]]^sup ~" mean, wIQR " ~ .(actualQQ[1]) ~ "-" ~ .(actualQQ[2]) ),xlab= xlab_i, ylab = ylab_i,...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})

      legend_info<-hist(zzBackgroud_max,plot=FALSE,breaks=20)
      which_col_leg <- which(legend_info$counts>0)
      par(mar=c(5,0.2,5,2.8))
      image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
      axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))

      par(mar=c(4,4,4,0.2) + 0.2)
      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=bkg_col,
            main=bquote(P[Psi[.(i)]]^inf ~" mean, wIQR " ~ .(actualQQ[1]) ~ "-" ~ .(actualQQ[2]) ),xlab= xlab_i, ylab = ylab_i,...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})

      legend_info<-hist(zzBackgroud_min,plot=FALSE,breaks=20)
      which_col_leg <- which(legend_info$counts>0)
      par(mar=c(5,0.2,5,2.8))
      image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
      axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))
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


      for(i in seq(num_Psi)){

#        whichVars<-as.logical(apply(allPsi[[i]]==1,2,sum))
#        pp <- bivProf$kmModel@X[,whichVars]

        wwVV<- apply(allPsi[[i]],2,function(x)return(x>0))

        pp<-bivProf$kmModel@X%*%t(allPsi[[i]])

        xlab_i <- paste(allPsi[[i]][1,wwVV[1,]],"*",plot_options$coord_names[wwVV[1,]],collapse = " + ")
        ylab_i <- paste(allPsi[[i]][2,wwVV[2,]],"*",plot_options$coord_names[wwVV[2,]],collapse = " + ")

        bck_name <- "mean"
        if(trueEvals<3){
          zzBackgroud_max<-matrix(allRes$res$max[[i]],ncol=dd_eta)
          zzBackgroud_min<-matrix(allRes$res$min[[i]],ncol=dd_eta)
          bck_name <- "mean"
        }else{
          propPoints<-getPointProportion(pp=pp,xBins = Design[1:(dd_eta),i],
                                         yBins = Design[(dd_eta+1):(2*dd_eta),i],whichAbove = whichAbove,plt = FALSE)

          zzBackgroud_max<-propPoints$freq
          zzBackgroud_min<-propPoints$freq
          bck_name <- "frequency"
        }
#        layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
        layout(matrix(1:4,nrow=1),widths=c(0.45,0.05,0.45,0.05))
        par(mar=c(4,4,4,0.2) + 0.2)
        image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=bkg_col,
              main=bquote(P[Psi[.(i)]]^sup ~ .(names(bivProf$res_UQ$bound$bound)[j]) ~ " bound," ~ " "~ .(bck_name)),xlab= xlab_i, ylab = ylab_i,...)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
        switch (as.character(trueEvals),
                "0" = {},
                "1" = {points(pp,pch=3)},
                {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})

        legend_info<-hist(zzBackgroud_max,plot=FALSE,breaks=20)
        which_col_leg <- which(legend_info$counts>0)
        par(mar=c(5,0.2,5,2.8))
        image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
        axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))

        par(mar=c(4,4,4,0.2) + 0.2)
        image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=bkg_col,
              main=bquote(P[Psi[.(i)]]^inf ~ .(names(bivProf$res_UQ$bound$bound)[j]) ~ " bound," ~ " "~ .(bck_name)),xlab= xlab_i, ylab = ylab_i,...)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T)
        contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=1.4)
        switch (as.character(trueEvals),
                "0" = {},
                "1" = {points(pp,pch=3)},
                {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})

        legend_info<-hist(zzBackgroud_min,plot=FALSE,breaks=20)
        which_col_leg <- which(legend_info$counts>0)
        par(mar=c(5,0.2,5,2.8))
        image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
        axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))
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

    for(i in seq(num_Psi)){

#      whichVars<-as.logical(apply(allPsi[[i]]==1,2,sum))
#      pp <- bivProf$kmModel@X[,whichVars]

      wwVV<- apply(allPsi[[i]],2,function(x)return(x>0))

      pp<-bivProf$kmModel@X%*%t(allPsi[[i]])

      xlab_i <- paste(allPsi[[i]][1,wwVV[1,]],"*",plot_options$coord_names[wwVV[1,]],collapse = " + ")
      ylab_i <- paste(allPsi[[i]][2,wwVV[2,]],"*",plot_options$coord_names[wwVV[2,]],collapse = " + ")

      zzBackgroud_max<-matrix(ddRange$max[,i],ncol=dd_eta)
      zzBackgroud_min<-matrix(ddRange$min[,i],ncol=dd_eta)

#      layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
      layout(matrix(1:4,nrow=1),widths=c(0.45,0.05,0.45,0.05))
      par(mar=c(4,4,4,0.2) + 0.2)
      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=bkg_col,
            main=bquote(P[Psi[.(i)]]^sup ~ "mean, bound range"),xlab= xlab_i, ylab = ylab_i,...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})


      legend_info<-hist(zzBackgroud_max,plot=FALSE,breaks=20)
      which_col_leg <- which(legend_info$counts>0)
      par(mar=c(5,0.2,5,2.8))
      image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
      axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))

      par(mar=c(4,4,4,0.2) + 0.2)
      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=bkg_col,
            main=bquote(P[Psi[.(i)]]^inf ~ "mean, bound range"),xlab= xlab_i, ylab = ylab_i,...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})

      legend_info<-hist(zzBackgroud_min,plot=FALSE,breaks=20)
      which_col_leg <- which(legend_info$counts>0)
      par(mar=c(5,0.2,5,2.8))
      image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
      axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))
    }

    if(plot_options$save)
      dev.off()


    ### Plot the bivariate profile mean as contour and sigma_T^Delta as background
    if(plot_options$save)
      pdf(file = paste(plot_options$folderPlots,"profUQ_sigmaD",plot_options$id_save,".pdf",sep=""),width = 18,height = 9,pointsize = 15)

    for(i in seq(num_Psi)){

#      whichVars<-as.logical(apply(allPsi[[i]]==1,2,sum))
#      pp <- bivProf$kmModel@X[,whichVars]

      wwVV<- apply(allPsi[[i]],2,function(x)return(x>0))

      pp<-bivProf$kmModel@X%*%t(allPsi[[i]])

      xlab_i <- paste(allPsi[[i]][1,wwVV[1,]],"*",plot_options$coord_names[wwVV[1,]],collapse = " + ")
      ylab_i <- paste(allPsi[[i]][2,wwVV[2,]],"*",plot_options$coord_names[wwVV[2,]],collapse = " + ")

      zzBackgroud_max<-matrix(sqrt(pmax(bivProf$res_UQ$bound$mean_var_D$var$res$max[,i],0)),ncol=dd_eta)
      zzBackgroud_min<-matrix(sqrt(pmax(bivProf$res_UQ$bound$mean_var_D$var$res$min[,i],0)),ncol=dd_eta)


#      zzBackgroud_max<-matrix(ddRange$max[,i],ncol=dd_eta)
#      zzBackgroud_min<-matrix(ddRange$min[,i],ncol=dd_eta)

#      layout(matrix(1:2,nrow=1),widths=c(0.5,0.5))
      layout(matrix(1:4,nrow=1),widths=c(0.45,0.05,0.45,0.05))
      par(mar=c(4,4,4,0.2) + 0.2)
      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_max,col=bkg_col,
            main=bquote(P[Psi[.(i)]]^sup ~ "mean, " ~(sigma[T]^tilde(Delta)) ),xlab= xlab_i, ylab = ylab_i,...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$max[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})


      legend_info<-hist(zzBackgroud_max,plot=FALSE,breaks=20)
      which_col_leg <- which(legend_info$counts>0)
      par(mar=c(5,0.2,5,2.8))
      image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
      axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))

      par(mar=c(4,4,4,0.2) + 0.2)
      image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackgroud_min,col=bkg_col,
            main=bquote(P[Psi[.(i)]]^inf ~ "mean, " ~(sigma[T]^tilde(Delta)) ),xlab= xlab_i, ylab = ylab_i,...)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),nlevels = 10,add=T,lwd = 1.5,labcex=1)
      contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(allRes$res$min[[i]],ncol=dd_eta),levels = threshold,add=T,col=2,lwd=2,labcex=1,lty=2)
      switch (as.character(trueEvals),
              "0" = {},
              "1" = {points(pp,pch=3)},
              {points(pp[-whichAbove,],pch=3,col="blue3");points(pp[whichAbove,],pch=3,col=3)})

      legend_info<-hist(zzBackgroud_min,plot=FALSE,breaks=20)
      which_col_leg <- which(legend_info$counts>0)
      par(mar=c(5,0.2,5,2.8))
      image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
      axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))
    }

    if(plot_options$save)
      dev.off()

  }


}
