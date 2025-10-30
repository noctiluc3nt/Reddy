#' Multiresolution Decomposition (MRD) according to Vickers and Mahrt, 2003
#'
#'@description Calculates multiresolution decomposition (MRD) according to Vickers and Mahrt, 2003
#'@param var1 timeseries of a variable
#'@param var2 timeseries of another variable to calculate the cospectrum of \code{var1} and \code{var2}, optional (default is \code{NULL})
#'@param time_res time resolution of the given timeseries in seconds (e.g., \code{time_res = 0.05} for 20 Hz)
#'@param plot logical, should the MRD spectrum be plotted? default \code{plot=TRUE}
#'@param ... arguments passed to plot function
#'@return MRD in form of data.frame with columns: index, m, scale, time, mean, median, q25, q75
#'
#'@export
#'
#'@examples
#'series=c(1,3,2,5,1,2,1,3) #example used in Vickers and Mahrt, 2003
#'calc_mrd(series)
#'
calc_mrd = function(var1,var2=NULL,time_res=0.05,plot=TRUE,...) {
    nf=length(var1)
    #calc exponent M
    pot=2^(0:100)
    M=log(pot[sum(pot<=nf)],2)
    #interpolate on new length 2^M
    var1=approx(var1,n=2^M)$y    
    if (!is.null(var2)) var2=approx(var2,n=2^M)$y
    #just reducing size of vector
    #var1=var1[1:2^M]
    #if (!is.null(var2)) var2=var2[1:2^M]
    #allocate output
    mrd_mean=array(NA,dim=M+1)
    mrd_median=array(NA,dim=M+1)
    mrd_q25=array(NA,dim=M+1)
    mrd_q75=array(NA,dim=M+1)
    mrd_sd=array(NA,dim=M+1)
    #mrd
    ms=M:0
    seg_mean=as.numeric(lapply(var1,mean,na.rm=TRUE))
    var1=var1-mean(var1,na.rm=TRUE)
    if (!is.null(var2)) var2=var2-mean(var2,na.rm=TRUE)
    for (i in 0:M) {
        if (is.null(var2)) {
            var1=matrix(var1,ncol=2^(M-i),byrow=TRUE)
            wn = apply(var1,1,mean,na.rm=TRUE)
            mrd_mean[M-i+1]=mean(wn^2,na.rm=TRUE)
            mrd_median[M-i+1]=median(wn^2,na.rm=TRUE)
            mrd_q25[M-i+1]=quantile(wn^2,0.25,na.rm=TRUE)
            mrd_q75[M-i+1]=quantile(wn^2,0.75,na.rm=TRUE)
            mrd_sd[M-i+1]=sd(wn^2,na.rm=T)
            var1=c(t(var1-wn))
        } else {
            var1=matrix(var1,ncol=2^(M-i),byrow=TRUE)
            wn = apply(var1,1,mean,na.rm=TRUE)
            var2=matrix(var2,ncol=2^(M-i),byrow=TRUE)
            wn2 = apply(var2,1,mean,na.rm=TRUE)
            mrd_mean[M-i+1]=mean(wn*wn2,na.rm=TRUE)
            mrd_median[M-i+1]=median(wn*wn2,na.rm=TRUE)
            mrd_q25[M-i+1]=quantile(wn*wn2,0.25,na.rm=TRUE)
            mrd_q75[M-i+1]=quantile(wn*wn2,0.75,na.rm=TRUE)
            mrd_sd[M-i+1]=sd(wn*wn2,na.rm=TRUE)
            var1=c(t(var1-wn))
            var2=c(t(var2-wn2))
        }
    }
    mrd_times=time_res*(2^(M:0))
    out=data.frame("index"=1:(M+1),"m"=M:0,"scale"=2^(M:0),"time"=mrd_times,"mean"=rev(mrd_mean),"median"=rev(mrd_median),"q25"=rev(mrd_q25),"q75"=rev(mrd_q75))
    if (plot==TRUE) plot_mrd(out, ...)
    return(out)
}

#' Suggestion of Averaging Time from Multiresolution Decomposition
#'
#'@description Suggested averaging time, i.e. first zero-crossing of MRD
#'@param mrd_out an object returned from \code{calc_mrd}
#'@return vector or scalar containing all zero crossing
#'@export
#'
#'@importFrom rootSolve uniroot.all
#'
#'@examples
#'set.seed(5)
#'series=rnorm(2^10)
#'mrd_test=calc_mrd(c(series))
#'suggest_avgtime_from_mrd(mrd_test)
#'
suggest_avgtime_from_mrd=function(mrd_out) {
    zerox=rootSolve::uniroot.all(approxfun(mrd_out$time,mrd_out$median),interval=range(mrd_out$time))
    return(zerox)
}

#' Plotting Multiresolution Decomposition
#'
#'@description Plots multiresolution decomposition (MRD)
#'@param mrd_out an object returned from \code{calc_mrd}
#'@param suggest_avgtime logical, should the suggested averaging time be plotted? default \code{TRUE}
#'@param ... arguments passed to plot function
#'@return creates a plot of MRD with logarithmic time scale (no return)
#'@export
#'
#'@examples
#'set.seed(5)
#'series=rnorm(2^10)
#'mrd_test=calc_mrd(c(series))
#'plot_mrd(mrd_test)
#'
plot_mrd=function(mrd_out,suggest_avgtime=TRUE,...) {
    if (!exists("ylab")) { ylab="MRD" }
    if (!exists("xlab")) { xlab="time [s]"}
    if (!exists("ylim")) { ylim=c(min(mrd_out$q25[!is.na(mrd_out$q25)]),max(mrd_out$q75[!is.na(mrd_out$q75)])) }
    plot(log10(mrd_out$time),mrd_out$median,lwd=2,pch=4,col=4,xlab=xlab,ylab=ylab,ylim=ylim,xaxt="n",type="b",...)
    suppressWarnings({
	    #arrows(log10(mrd_out$time),mrd_out$median,log10(mrd_out$time),mrd_out$q75,lwd=1,col=4,angle=90,length=0.05)
	    #arrows(log10(mrd_out$time),mrd_out$median,log10(mrd_out$time),mrd_out$q25,lwd=1,col=4,angle=90,length=0.05)
        shade_between(log10(mrd_out$time),log10(mrd_out$time),mrd_out$q75,mrd_out$q25,col=rgb(0,0,0.9,0.1),lty=0)
    })
	points(log10(mrd_out$time),mrd_out$mean,lwd=2,col=1,type="b",pch=20)
	axis(1,at=-1:4,labels=10^(-1:4))
	legend("topleft",legend=c("mean","median"),col=c(1,4),lty=0,pch=c(20,4),lwd=2)
	abline(v=log10(60),lty=3)
	abline(v=log10(60*30),lty=3)
    abline(h=0,lty=3)
    if (suggest_avgtime==TRUE) {
        zerox=suggest_avgtime_from_mrd(mrd_out)[1]
        if (!is.na(zerox)) {
            abline(v=min(log10(zerox)),col="orangered")
            print(paste("suggested averaging time (i.e. first zero-crossing):",round(zerox/60,2),"min"))
        } else {
            print("No zero-crossing of the MRD (co)spectrum was detected.")
        }
    }
}
