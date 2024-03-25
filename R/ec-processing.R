#' Despiking
#'
#'@description Three despiking method based on 1) pre-defined thresholds, 2) median deviation (mad) test and 3) skewness and kurtosis
#'@param series timeseries that shall be despiked
#'@param thresholds vector with two elements representing lower and upper bounds for despiking (pre-defined thresholds), 'NA' means that the respective bound is not used
#'@param mad_factor factor for the mad test, default 'mad_factor = 10'
#'@param threshold_skewness threshold for skewness test, default 'threshold_skewness = 2'
#'@param threshold_kurtosis threshold for kurtosis test, default 'threshold_kurtosis = 8'
#'
#'@return despiked timeseries
#'@export
#'
#'@importFrom pracma detrend
#'
#'@examples
#'set.seed(5)
#'ts1=rnorm(100)
#'despiking(ts1,thresholds=c(-1,1))
#'
#'ts2=rexp(1000)
#'despiking(ts2)
#'
despiking = function(series,thresholds=c(NA,NA),mad_factor=10,threshold_skewness=2,threshold_kurtosis=8) {
	#despiking based on predefined limits
    if (sum(is.na(thresholds))==0) {
        pass = (series>thresholds[1] & series<thresholds[2])
	    series[!pass] = NA
    } else if (is.na(thresholds[1])) {
        pass = (series<thresholds[2])
	    series[!pass] = NA
    } else if (is.na(thresholds[2])) {
        pass = (series>thresholds[1])
	    series[!pass] = NA
    }
    #despiking based on median deviation test
    med=median(series,na.rm=TRUE)
	mad=median(abs(series-med),na.rm=TRUE)
	pass=(abs(series-med) <= mad_factor*mad) #pass criterion
	series[!pass] = NA
    #despiking based on skewness and kurtosis
    seriesLDT=pracma::detrend(series,tt="linear") #linear detrending to eliminate trends (departures from stationarity) -> would influnece higher moments
	skewness=mean(seriesLDT^3)/sd(seriesLDT)^3
	kurtosis=mean(seriesLDT^4)/sd(seriesLDT)^4
	pass=(abs(skewness)<threshold_skewness & kurtosis<threshold_kurtosis)
	series[!pass] = NA
    return(series)
}


#' Double rotation
#'
#'@description Double rotation
#'@param u u-wind (levelled sonic)
#'@param v v-wind (levelled sonic)
#'@param w w-wind (levelled sonic)
#'
#'@return list containing the wind in a natural coordinate system (streamwise, crosswise, vertical) and the two rotation angles theta and phi
#'@export
#'
#'@examples
#'
#'
rotate_double = function(u,v,w) {
	#horizontal
	theta=atan2(mean(v,na.rm=T),mean(u,na.rm=T))
	u1=u*cos(theta) + v*sin(theta)
    v1=-u*sin(theta) + v*cos(theta)
    w1=w
	#vertical
	phi=atan2(mean(w1,na.rm=T),mean(u1,na.rm=T))
	u2=u1*cos(phi) + w1*sin(phi)
    v2=v1
    w2=-u1*sin(phi)+w1*cos(phi)
	return(list("u"=u2,"v"=v2,"w"=w2,"theta"=theta*180/pi,"phi"=phi*180/pi))
}


#' Planar fit rotation
#'
#'@description Planar fit rotation
#'@param u u-wind (levelled sonic)
#'@param v v-wind (levelled sonic)
#'@param w w-wind (levelled sonic)
#'
#'@return 
#'@export
#'
#'@examples
#'
#'
rotate_planar = function(u,v,w) {
	#TODO
	theta=atan2(mean(v),mean(u))
	rot1=matrix(c(cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1), nrow=3,ncol=3,byrow=TRUE) #B bzw. M^T
	wind1=c(u,v,w)%*%rot1
	return(list("wind"=wind1,"theta"=theta*180/pi))
}


#' Stationarity Flag
#'
#'@description Stationarity Flag according to Foken and Wichura, 1996 based on the assumption that the covariance of two variables ('var1' and 'var2', one usually representing vertical velocity) calculated for blocks (of length nsub) do not differ to much from the total covariance
#'@param var1 variable 1 
#'@param var2 variable 2 (same length as 'var1', usually either 'var1' or 'var2' represent vertical velocity)
#'@param nsub number of elements used for subsampling (nsub < length(var1)) 
#'
#'@return 
#'@export
#'
#'@examples
#'set.seed(5)
#'ts1=rnorm(30)
#'ts2=rnorm(30)
#'flag_stationarity(ts1,ts2)
#'
flag_stationarity = function(var1,var2,nsub=6) {
    if (length(var1) != length(var2)) {
        print("ERROR: var1 and var2 have to be of equal length.")
    }
    nint=length(var1)%/%nsub
    if (nint<=1) {
        print("WARNING: nsub is chosen to large.")
    }
    rnfs=array(NA,dim=nint)
    cov_complete=cov(var1,var2,use="pairwise.complete.obs")
    for (i in 1:nint) {
        isub=((i-1)*nsub+1):(i*nsub)
        cov_sub=cov(var1[isub],var2[isub],use="pairwise.complete.obs")
        rnfs[i]=abs((cov_complete-cov_sub)/cov_complete)
    }
    rnf=mean(rnfs,na.rm=T)
    flag=ifelse(rnf<0.3,0,ifelse(rnf<1,1,2))
    return(flag)
}


#' Vertical Velocity Flag
#'
#'@description Vertical Velocity Flag according to Mauder et al., 2013: After (planar fit) rotation the vertical velocity should vanish, this flag flags high remaining vertical velocities.
#'@param w vertical velocity 
#'
#'@return 
#'@export
#'
#'@examples
#'
flag_w = function(w) {
    w=abs(w)
    flag=ifelse(w<0.1,0,ifelse(1<0.15,1,2))
    return(flag)
}

#' Flow Distortion Flag and Wind Constancy Ratio
#'
#'@description Flow Distortion Flag according to Mauder et al., 2013: Wind coming from (pre-defined) directions blocked by the measurement device is flaged with 2 (for wind speeds greater than 0.1 assuming that during calm wind the wind direction is not well-defined). The wind constancy ratio is calculated to quantify the variability of horizontal wind direction according to Mahrt, 1999.
#'@param u u-wind (levelled sonic)
#'@param v v-wind (levelled sonic)
#'@param dir_blocked vector containing the lower and upper bound of the blocked wind sector in degrees (e.g., 'dir_blocked=c(30,60)')
#'@param threshold_cr threshold for constancy ratio (default 'threshold_cr=0.9', may be adapted to used data set)
#'
#'@return 
#'@export
#'
#'@examples
#'
flag_distortion = function(u,v,dir_blocked=c(30,60),threshold_cr=0.9) {
    if (length(u) != length(v)) {
        print("ERROR: u and v have to be of equal length.")
    }
    #horizontal wind speed
    ws=sqrt(mean(u,na.rm=T)^2+mean(v,na.rm=T)^2)
    #constancy ratio cr
    cr=sqrt(sum(u^2,na.rm=T)+sum(v^2,na.rm=T))/(ws)
    #flow distortion flag considering cr
    if (!is.na(ws) & !is.na(cr)) {
        if (ws>0.1 & cr>threshold_cr) {
            wd=atan2(mean(v),mean(u))
            flag=ifelse(wd>=dir_blocked[1] & wd<=dir_blocked[2],2,0)
        } else {
            flag=0
        }
    } else {
        flag=2
    }   
    return(flag)
}

#' Integral Turbulence Characteristics Flag 
#'
#'@description Integral Turbulence Characteristics Flag: Tests the consistency with Monin-Obukhov similarity theory using the scaling functions from Panofsky and Dutton, 1984
#'@param sigma_w standard deviation of vertical velocity
#'@param ustar friction velocity
#'@param zeta stability parameter zeta = z/L
#'
#'@return 
#'@export
#'
#'@examples
#'
flag_most = function(sigma_w,ustar,zeta) {
    parameterized=1.3*(1+2*abs(zeta))^(1/3) #sigma_w/ustar parametrized according to scaling function based on zeta
    itc=abs((sigma_w/ustar-parameterized)/parameterized)
    flag=ifelse(itc<0.3,0,ifelse(itc<0.8,1,2))
    return(flag)
}

#' SND correction
#'
#'@description SND correction of sensible heat flux
#'@param u u-wind (levelled sonic)
#'@param v v-wind (levelled sonic)
#'@param w w-wind (levelled sonic)
#'@param Ts temperature (sonic temperature or corrected temperature)
#'@param q specific humidity (if measured, default NULL)
#'@param A constant used in SND correction, default 'A = 7/8'
#'@param B constant used in SND correction, default 'B = 7/8'
#'
#'@return 
#'@export
#'
#'@examples
#'
SNDcorrection = function(u,v,w,Ts,q=NULL,A=7/8,B=7/8) {
    #calculation of respective covariances
	not_na=!is.na(w)&!is.na(Ts)
	covar_wTs = cov(w[not_na],Ts[not_na]) 
	not_na=!is.na(w)&!is.na(u)
	covar_uw = cov(u[not_na],w[not_na]) 
	not_na=!is.na(w)&!is.na(v)
	covar_vw = cov(v[not_na],w[not_na])
	ubar=mean(u,na.rm=TRUE)
	vbar=mean(v,na.rm=TRUE)
	Tsbar=mean(Ts,na.rm=TRUE)
	if (!is.null(q)) { #considering q
		not_na=!is.na(w)&!is.na(q)
		covar_qw=covar(q[not_na],w[not_na])
		return(covar_wTs - 0.51*covar_qw + 2*Tsbar/clight()^2*(A*ubar*covar_uw + B*vbar*covar_vw))
	}
    #without q
	return(covar_wTs + 2*Tsbar/clight()^2*(A*ubar*covar_uw + B*vbar*covar_vw))
}

