#' discrete binning
#'
#'@description discrete binning of a variable var1 based on another variable var2 (e.g., the stability parameter zeta)
#'@param var1 vector, variable that should be binned
#'@param var2 vector, variable used for the binning
#'@param bins vector, providing the intervals of the bins of var2
#'@return matrix of dimension (length(bins)-1,4) with columns representing mean, median, 25%-quantile, 74%-quantile
#'@export
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
			    #out[i-1,3]=median(sub,na.rm=T)+sd(sub,na.rm=T)/sqrt(nm) #using the estimator uncertainty
			    #out[i-1,4]=median(sub,na.rm=T)-sd(sub,na.rm=T)/sqrt(nm)
			    out[i-1,3]=quantile(sub,0.25,na.rm=T)
			    out[i-1,4]=quantile(sub,0.75,na.rm=T)
	    	}
	    }
	    return(out)
    }
}