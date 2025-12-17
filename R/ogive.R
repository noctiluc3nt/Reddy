#' Calculates Ogive (cumulative distribution function from cospectrum) based on MRD
#'
#'@description Calculates Ogive from MRD spectram
#'@param mrd an object returned from \code{calc_mrd}
#'@param plot logical, should the ogive be plotted? default \code{plot=TRUE}
#'@param ... arguments passed to plot function
#'@return ogive as list containing index, m, scale, time, ogive based on mean, and ogive based on median with quantiles q25, q75
#'@export
#'
#'@examples
#'set.seed(5)
#'series=rnorm(2^10)
#'mrd_test=calc_mrd(c(series))
#'ogive_test=calc_ogive(mrd_test,plot=FALSE)
#'
calc_ogive=function(mrd,plot=TRUE,...) {
    nt=length(mrd$index)
    ogive_mean=array(NA,dim=c(nt))
    ogive_median=array(NA,dim=c(nt))
    ogive_q25=array(NA,dim=c(nt))
    ogive_q75=array(NA,dim=c(nt))
    for (i in 1:nt) {
        inds=(i:nt)
        ogive_mean[i]=sum(mrd$mean[inds])
        ogive_median[i]=sum(mrd$median[inds])
        ogive_q25[i]=sum(mrd$q25[inds])
        ogive_q75[i]=sum(mrd$q75[inds])
    }
    ogive=data.frame("index"=mrd$index,"m"=mrd$m,"scale"=mrd$scale,"time"=mrd$time,
        "mean"=ogive_mean,"median"=ogive_median,"q25"=ogive_q25,"q75"=ogive_q75)
    if (plot==TRUE) {
        tryCatch({
            plot_ogive(ogive,...)
        }, warning=function(e){
            message("An error occurred when plotting the ogive.")
        })
    }
    return(ogive)
}


#' Plotting Ogive
#'
#'@description Plots ogive
#'@param ogive an object returned from \code{calc_ogive}
#'@param ... arguments passed to plot function
#'@return creates a plot of an ogive with logarithmic time scale (no return)
#'@export
#'
#'@examples
#'set.seed(5)
#'series=rnorm(2^10)
#'mrd_test=calc_mrd(c(series))
#'ogive_test=calc_ogive(mrd_test)
#'plot_ogive(ogive_test)
#'
plot_ogive=function(ogive,...) {
    #that's in principle the same function as plot_mrd()
    if (!exists("ylab")) { ylab="Ogive" }
    if (!exists("xlab")) { xlab="time [s]"}
    if (!exists("ylim")) { ylim=c(min(ogive$q25[!is.na(ogive$q25)]),max(ogive$q75[!is.na(ogive$q75)])) }
    plot(log10(ogive$time),ogive$median,lwd=2,pch=4,col=4,xlab=xlab,ylab=ylab,ylim=ylim,xaxt="n",type="o",...)
    suppressWarnings({
	    #arrows(log10(ogive$time),ogive$median,log10(ogive$time),ogive$q75,lwd=1,col=4,angle=90,length=0.05)
	    #arrows(log10(ogive$time),ogive$median,log10(ogive$time),ogive$q25,lwd=1,col=4,angle=90,length=0.05)
        shade_between(log10(ogive$time),log10(ogive$time),ogive$q75,ogive$q25,col=rgb(0,0,0.9,0.1),lty=0)
    })
	points(log10(ogive$time),ogive$mean,lwd=2,col=1,type="o",pch=20)
	axis(1,at=-1:4,labels=10^(-1:4))
	legend("topleft",legend=c("mean","median"),col=c(1,4),lty=0,pch=c(20,4),lwd=2)
	abline(v=log10(60),lty=3)
	abline(v=log10(60*30),lty=3)
    abline(h=0,lty=3)
}