#' discrete binning
#'
#'@description discrete binning of a variable \code{var1} based on another variable \code{var2} (e.g., the stability parameter \code{zeta})
#'@param var1 vector, variable that should be binned
#'@param var2 vector, variable used for the binning
#'@param bins vector, providing the intervals of the bins of \code{var2}
#'@return matrix of dimension (\code{length(bins)-1,4}) with columns representing mean, median, q25, q75
#'@export
#'
#'@examples
#'zeta_bins=c(-10^(3:-3),10^(-3:3))
#'zeta_vals=rnorm(1000)
#'vals=runif(1000)
#'binned=binning(vals,zeta_vals,zeta_bins)
#'
binning=function(var1,var2,bins) {
    if (length(var1)!=length(var2)) {
        print("ERROR: var1 and var2 have to have the same length")
    } else {
        nbins=length(bins)
	    out=array(NA,dim=c(nbins-1,4))
	    for (i in 2:nbins) {
		    sub=var1[var2>bins[i-1] & var2<bins[i]]
		    if (any(!is.na(sub))) {
			    out[i-1,1]=mean(sub,na.rm=T)
			    out[i-1,2]=median(sub,na.rm=T)
			    #nm=length(sub)
			    #out[i-1,3]=median(sub,na.rm=T)+sd(sub,na.rm=T)/sqrt(nm) #using the estimation uncertainty
			    #out[i-1,4]=median(sub,na.rm=T)-sd(sub,na.rm=T)/sqrt(nm)
			    out[i-1,3]=quantile(sub,0.25,na.rm=T)
			    out[i-1,4]=quantile(sub,0.75,na.rm=T)
	    	}
	    }
	    return(out)
    }
}

#' Shifting two timeseries to match maximum cross-correlation
#'
#'@description Shifts two timeseries to match their maximum cross-correlation
#'@param var1 vector, first timeseries
#'@param var2 vector, second timeseries
#'@param plot logical, should the cross-correlation be plotted? default \code{plot = TRUE}
#'@return a matrix cotaining timeseries \code{var1} and \code{var2} as columns after shifting to the maximum cross-correlation
#'@export
#'
#'@examples
#'ts1=runif(10)
#'ts2=c(1,1,ts1)
#'shifted=shift2maxccf(ts1,ts2)
#'
shift2maxccf=function(var1,var2,plot=TRUE) {
	#equalize lengths of timeseries
	n1=length(var1)
	n2=length(var2)
	if (n1<n2) {
		var1=c(var1,rep(NA,n2-n1))
	} else if (n2<n1) {
		var2=c(var2,rep(NA,n1-n2))
	}
	n=max(n1,n2)
	#calc cross-correlation
	notna=(!is.na(var1) & !is.na(var2))
    cc=ccf(var1[notna],var2[notna])
	maxcc=max(cc$acf) #max positive cross-correlation
	lag=cc$lag[which(cc$acf==max(cc$acf))] #respective time lag
	if (plot == TRUE) {
		plot(cc)
		points(lag,maxcc,col=2,lwd=2,cex=2)
		print(paste("lag:",lag,"\n maximum cross-correlation:",maxcc))
	}
	#shift timeseries
	mat=array(NA,dim=c(n+abs(lag),2))
	if (lag>0) { #positive lag = var1 leads var2
		mat[,1]=c(var1,rep(NA,lag))
		mat[,2]=c(rep(NA,lag),var2[1:n])
	} else { #negative lag = var2 leads var1
		mat[,1]=c(rep(NA,abs(lag)),var1[1:n])
		mat[,2]=c(var2,rep(NA,abs(lag)))
	}
	return(mat)
}

#' deaccumulation 
#'
#'@description hourly deaccumulation, e.g. for fluxes from model output
#'@param dat vector (with dimension time) or array (with dimension x, y, time)
#'@param factor factor for unit and sign conversion, default: \code{factor = -1/3600} for converting hour to seconds and adapting the sign convention
#'@return vector or array hourly deaccumulated (same dimension as input)
#'@export
#'
deaccumulate1h=function(dat,factor=-1/3600) {
	if (is.vector(dat)) {
		n=length(dat)
		out=dat
    	for (i in 2:n) {
        	out[i]=dat[i]-dat[i-1] #deaccumulate hourly
    	}
	} else {
		dims=dim(dat)
		n=dims[3]
		field_out=dat
   		for (i in 2:n) {
        	out[,,i]=dat[,,i]-dat[,,i-1] #deaccumulate hourly
    	}
	}
    return(out*factor) #time unit conversion
}


#' gap-filling
#'
#'@description gap-filling of a timeseries based on linear or constant interpolation
#'@param var timeseries, where NA indicates missing values that should be filled
#'@param method interpolation method, can be either \code{method = "linear"} for linear interpolation (default) or \code{method = "constant"} for constant interpolation
#'@param nmissing number of allowed missing values, default \code{nmissing = 4}
#'@return gap-filled timeseries
#'@export
#'
#'@examples
#'ts1=c(1,2,NA,0)
#'gapfilling(ts1) #1,2,1,0
#'gapfilling(ts1,method="constant") #1,2,2,0
#'gapfilling(ts1,nmissing=0) #too many missing values
#'
gapfilling=function(var,nmissing=4,method="linear") {
    n=length(var)
	nm=sum(is.na(var))
	if (nm <= nmissing) {
		var_interpolated = approx(1:n,var,xout=1:n,method=method)
		return(var_interpolated$y)
	} else {
		warning("The timeseries contains more missing values than desired (specified in nmissing).")
	}
}


#' accumulating / averaging
#'
#'@description averaging of a timeseries
#'@param var timeseries
#'@param tres1 time resolution [s] of the given timeseries \code{var}, default \code{tres1 = 0.05} (for 20 Hz)
#'@param tres2 desired time resolution(s) [s] of the averaged timeseries (scalar or vector), default \code{tres2 = c(1,10,30)*60} (for 1, 10 and 30 minutes)
#'@return list containing mean and standard deviation of the timeseries for the desired time interval(s)
#'@export
#'@importFrom RcppRoll roll_mean roll_sd
#'
#'@examples
#'ts=rnorm(30*60*20) #30 minutes of 20 Hz measurements
#'averaging(ts)
#'
averaging=function(var,tres1=0.05,tres2=c(1,10,30)*60) {
	n=length(var)
	nt=length(tres2)
    maxt=max(tres2)
	if (n<(maxt/tres1)) warning("The timeseries is not long enough for the desired (maximum) averaging time.")
	nav=tres2/tres1 #number of values to be averaged
	averaged_mean=list()
	averaged_sd=list()
	for (i in 1:nt) {
		vari=RcppRoll::roll_mean(var,nav[i])[seq(1,n,nav[i])]
		averaged_mean[[i]]=vari
		vari=RcppRoll::roll_sd(var,nav[i])[seq(1,n,nav[i])]
		averaged_sd[[i]]=vari
	}
	averaged=list("mean"=averaged_mean,"sd"=averaged_sd)
	averaged$averaging_time_min=tres2/60
	averaged$number_of_averaged_values=nav
	return(averaged)
}