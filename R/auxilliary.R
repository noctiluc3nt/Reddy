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
