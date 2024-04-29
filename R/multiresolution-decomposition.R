#' Multiresolution Decomposition (MRD) according to Vickers and Mahrt, 2003
#'
#'@description Calculates multiresolution decomposition (MRD) according to Vickers and Mahrt, 2003
#'@param var1 timeseries of a variable
#'@param var2 timeseries of another variable to calculate cospectrum of \code{var1} and \code{var2}, optional (default is \code{NULL})
#'@param time_res time resolution of the given timeseries in seconds (e.g., \code{time_res = 0.05} for 20 Hz)
#'@return MRD in form of data frame containing the columns: index, scale, time, mean, median, q25, q75
#'@export
#'
#'@examples
#'series=c(1,3,2,5,1,2,1,3) #example used in Vickers and Mahrt, 2003
#'calc_mrd(series)
#'
calc_mrd = function(var1,var2=NULL,time_res=0.05) {
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
    seg_mean=as.numeric(lapply(var1,mean,na.rm=T))
    var1=var1-mean(var1,na.rm=T)
    if (!is.null(var2)) var2=var2-mean(var2,na.rm=T)
    for (i in 0:M) {
        if (is.null(var2)) {
            var1=matrix(var1,ncol=2^(M-i),byrow=T)
            wn = apply(var1,1,mean,na.rm=T)
            mrd_mean[M-i+1]=mean(wn^2,na.rm=T)
            mrd_median[M-i+1]=median(wn^2,na.rm=T)
            mrd_q25[M-i+1]=quantile(wn^2,0.25,na.rm=T)
            mrd_q75[M-i+1]=quantile(wn^2,0.75,na.rm=T)
            mrd_sd[M-i+1]=sd(wn^2,na.rm=T)
            var1=c(t(var1-wn))
        } else {
            var1=matrix(var1,ncol=2^(M-i),byrow=T)
            wn = apply(var1,1,mean,na.rm=T)
            var2=matrix(var2,ncol=2^(M-i),byrow=T)
            wn2 = apply(var2,1,mean,na.rm=T)
            mrd_mean[M-i+1]=mean(wn*wn2,na.rm=T)
            mrd_median[M-i+1]=median(wn*wn2,na.rm=T)
            mrd_q25[M-i+1]=quantile(wn*wn2,0.25,na.rm=T)
            mrd_q75[M-i+1]=quantile(wn*wn2,0.75,na.rm=T)
            mrd_sd[M-i+1]=sd(wn*wn2,na.rm=T)
            var1=c(t(var1-wn))
            var2=c(t(var2-wn2))
        }
    }
    out=data.frame("index"=1:(M+1),"m"=M:0,"scale"=2^(M:0),"time"=time_res*(2^(M:0)),"mean"=mrd_mean,"median"=mrd_median,"q25"=mrd_q25,"q75"=mrd_q75)
    return(out)
}


#' Plotting Multiresolution Decomposition
#'
#'@description Plots multiresolution decomposition (MRD)
#'@param mrd_out an output object from \code{calc_mrd}
#'@param ... parameters passed to plot function
#'@return creates a plot of MRD with logarithmic time scale (no return)
#'@export
#'
#'@examples
#'set.seed(5)
#'series=rnorm(2^10)
#'mrd_test=calc_mrd(c(series))
#'plot_mrd(mrd_test)
#'
plot_mrd=function(mrd_out,...) {
    if (!exists("ylab")) { ylab="MRD" }
    if (!exists("xlab")) { xlab="averaging time [s]"}
    if (!exists("ylim")) { ylim=c(min(mrd_out$q25[!is.na(mrd_out$q25)]),max(mrd_out$q75[!is.na(mrd_out$q75)])) }
    plot(log10(mrd_out$time),mrd_out$median,lwd=2,pch=3,col=4,xlab=xlab,ylab=ylab,ylim=ylim,xaxt="n",type="b",...)
    suppressWarnings({
	    arrows(log10(mrd_out$time),mrd_out$median,log10(mrd_out$time),mrd_out$q75,lwd=1,col=4,angle=90,length=0.05)
	    arrows(log10(mrd_out$time),mrd_out$median,log10(mrd_out$time),mrd_out$q25,lwd=1,col=4,angle=90,length=0.05)
    })
	points(log10(mrd_out$time),mrd_out$mean,lwd=2,col=1,type="b")
	axis(1,at=-1:4,labels=10^(-1:4))
	legend("topleft",legend=c("mean","median"),col=c(1,4),lty=0,pch=c(1,3),lwd=2)
	abline(v=log10(60),lty=3)
	abline(v=log10(60*30),lty=3)
    abline(h=0,lty=3)
}