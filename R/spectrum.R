#' Turbulence Spectrum of Timeseries
#'
#'@description Calculates and plots the averaged turbulence spectrum (as wrapper of rbase::spectrum)
#'@param ts timeseries
#'@param nbins number of bins used to average the spectrum, default \code{nbins=100}
#'@param plot should the spectrum be plotted? default \code{plot=TRUE}
#'
#'@return binned spectrum
#'@export
#'
calc_spectrum = function(ts,nbins=100,plot=TRUE) {
	s=spectrum(ts,plot=FALSE)
    bins=seq(log(min(s$freq,na.rm=TRUE)),log(max(s$freq,na.rm=T)),length.out=nbins)
    sbin=binning(s$spec,s$freq,10^bins)
    if (plot==TRUE) {
        plot(bins[2:nbins],log(sbin[,2]),pch=20,main="Spectrum",xlab="frequency",ylab="spectrum",xlim=range(bins[2:nbins][!is.na(sbin[,2])],na.rm=T),xaxt="n",yaxt="n")
        axis(1,at=(-10):0,10^(-10:0))
        axis(2,at=(-10):10,10^(-10:10))
        fit=lm(log(sbin[,2]) ~ bins[2:nbins])
        print(summary(fit))
        par(new=T)
        plot(bins,bins*(-5/3),type="l",col=4,lwd=2,lty=2,xlab="",ylab="",xaxt="n",yaxt="n",ylim=range(log(sbin[,2]),na.rm=T))
        legend("bottomleft",legend="-5/3 slope",col=4,lwd=2,lty=2)
    }
    return(cbind(bins[2:nbins],sbin[,2]))
}
