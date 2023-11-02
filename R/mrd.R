#' Multiresolution Decomposition (MRD) according to Vickers and Mahrt, 2003
#'
#'@description Calculates Multiresolution Decomposition (MRD) according to Vickers and Mahrt, 2003
#'@param var1 timeseries of a variable
#'@param var2 timeseries of another variable to calculate cospectrum of var1 and var2, optional (default is NULL)
#'@param time_res time resolution of the given timeseries in seconds (e.g., 0.05 for 20 Hz)
#'@return MRD in form of data frame containing the columns: index, scale, time, mean, median, q25, q75
#'@export
#'
#'@examples
#'series=c(1,3,2,5,1,2,1,3)
#'mrd(series)
#'
mrd = function(var1,var2=NULL,time_res=0.05) {
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
            mrd_mean[i]=mean(wn^2,na.rm=T)
            mrd_median[i]=median(wn^2,na.rm=T)
            mrd_q25[i]=quantile(wn^2,0.25,na.rm=T)
            mrd_q75[i]=quantile(wn^2,0.75,na.rm=T)
            mrd_sd[i]=sd(wn^2,na.rm=T)
            var1=c(t(var1-wn))
        } else {
            var1=matrix(var1,ncol=2^(M-i),byrow=T)
            wn = apply(var1,1,mean,na.rm=T)
            var2=matrix(var2,ncol=2^(M-i),byrow=T)
            wn2 = apply(var2,1,mean,na.rm=T)
            mrd_mean[i]=mean(wn*wn2,na.rm=T)
            mrd_median[i]=median(wn*wn2,na.rm=T)
            mrd_q25[i]=quantile(wn*wn2,0.25,na.rm=T)
            mrd_q75[i]=quantile(wn*wn2,0.75,na.rm=T)
            mrd_sd[i]=sd(wn*wn2,na.rm=T)
            var1=c(t(var1-wn))
            var2=c(t(var2-wn2))
        }
    }
    out=data.frame("index"=1:(M+1),"scale"=2^(M:0),"time"=time_res*(2^(M:0)),"mean"=mrd_mean,"median"=mrd_median,"q25"=mrd_q25,"q75"=mrd_q75)
    return(out)
}

