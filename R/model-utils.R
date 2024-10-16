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


#' Converts hybrid (terrain-following) sigma levels to physical heights 
#'
#'@description Converts hybrid (terrain-following) sigma levels to physical heights
#'@param hybrid scalar or vector, hybrid sigma levels
#'@param Tv virtual temperature
#'@return hybrid levels converted to physical height [m]
#'@export
#'
sigma2height=function(hybrid,Tv=273.15) {
    return(-R()*Tv/g()*log(hybrid))
}