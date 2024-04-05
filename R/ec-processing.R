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
#'@param c a three-dimensional vector containing the offset of u-, v-, w-wind
#'
#'@return 
#'@export
#'
#'@examples
#'
rotate_planar = function(u,v,w,c=c(0,0,0)) {
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
#'flag_stationarity(ts1,ts2,nsub=6)
#'
flag_stationarity = function(var1,var2,nsub=3000) {
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

#' SND and cross-wind correction of sensible heat flux
#'
#'@description SND and cross-wind correction of sensible heat flux: converts the buoyancy flux cov(w,Ts) (based on sonic temperature Ts) to sensible heat flux
#'@param u u-wind [m/s] (levelled sonic)
#'@param v v-wind [m/s] (levelled sonic)
#'@param w w-wind [m/s] (levelled sonic)
#'@param Ts temperature [K] (sonic temperature or corrected temperature)
#'@param q specific humidity [kg/kg] (if measured by the sonic, default NULL)
#'@param A constant used in cross-wind correction, default 'A = 7/8' for CSAT3
#'@param B constant used in cross-wind correction, default 'B = 7/8' for CSAT3
#'
#'@return SND correction of sensible heat flux
#'@export
#'
SNDcorrection = function(u,v,w,Ts,q=NULL,A=7/8,B=7/8) {
    #calculation of respective covariances
	not_na=!is.na(w)&!is.na(Ts)
	cov_wTs = cov(w[not_na],Ts[not_na]) 
	not_na=!is.na(w)&!is.na(u)
	cov_uw = cov(u[not_na],w[not_na]) 
	not_na=!is.na(w)&!is.na(v)
	cov_vw = cov(v[not_na],w[not_na])
	ubar=mean(u,na.rm=TRUE)
	vbar=mean(v,na.rm=TRUE)
	Tsbar=mean(Ts,na.rm=TRUE)
	if (!is.null(q)) { #considering q
		not_na=!is.na(w)&!is.na(q)
		cov_qw=cov(q[not_na],w[not_na])
        #second term: SND correction, third term: cross-wind correction
		return(cov_wTs - 0.51*cov_qw + 2*Tsbar/clight()^2*(A*ubar*cov_uw + B*vbar*cov_vw))
	}
    #without q: only cross-wind correction
	return(cov_wTs + 2*Tsbar/clight()^2*(A*ubar*cov_uw + B*vbar*cov_vw))
}

#' WPL correction
#'
#'@description WPL correction: density correction for trace gas fluxes (i.e., converts volume- to mass-related quantity)
#'@param rho_w measured water vapor density [kg/m^3]
#'@param rho_c measured trace gas density [kg/m^3] (only if WPL-correction should be applied to another flux, e.g. CO2 flux, default NULL)
#'@param w w-wind [m/s] (levelled sonic)
#'@param Ts temperature [K] (sonic temperature or corrected temperature)
#'@param q specific humidity [kg/kg] (if measured, default NULL)
#'
#'@return WPL correction of respective flux
#'@export
#'
WPLcorrection = function(rho_w,rho_c=NULL,w,Ts,q) {
    #calculation of respective covariances
	not_na=!is.na(w)&!is.na(Ts)
	cov_wTs = cov(w[not_na],Ts[not_na]) 
	not_na=!is.na(w)&!is.na(rho_w)
	cov_wrhow = cov(w[not_na],rho_w[not_na]) 
    Ts_bar=mean(Ts,na.rm=T)
    q_bar=mean(q,na.rm=T)
    rho_w_bar=mean(rho_w,na.rm=T)
    if (is.null(rho_c)) { #water vapor flux
        return(1+1.61*q_bar)*(cov_wrhow+rho_w_bar/Ts_bar*cov_wTs) #with M_L/M_w = 1.61
    } else { #other trace gas flux
        not_na=!is.na(w)&!is.na(rho_c)
	    cov_wrhoc = cov(w[not_na],rho_c[not_na]) 
        rho_c_bar=mean(rho_c,na.rm=T)
        return(cov_wrhoc+1.61*rho_c_bar/rho_w_bar*cov_wrhow+(1+1.61*q_bar)*rho_c_bar/Ts_bar*cov_wTs)
    }
}


#' Unit conversion of "parts-per" to density (for closed-path gas analyzer)
#'
#'@description Unit conversion of "parts-per" to density (for closed-path gas analyzer)
#'@param ppt measurement in parts per thousand [ppt]
#'@param T_mean temperature [K]
#'@param pres pressure [Pa]
#'@param e water vapor pressure [Pa]
#'@param gas which gas? can be either "H2O", "CO2", "CH4" (if CO2/CH4 is selected, make sure that it's still in ppt and not ppm as usual)
#'
#'@return density of the gas [kg/m^3]
#'@export
#'
#'@examples
#'
ppt2rho = function(ppt,T_mean=288.15, pres = 101325, e = 0, gas="H2O") {
    Vd=Runiversal()*T_mean/(pres-e) #volume of dry air [m^3/mol]
    if (gas == "H2O") {
        return(ppt/1000*M_H2O()/Vd)
    } else if (gas == "CO2") {
        return(ppt/1000*M_CO2()/Vd)
    } else if (gas == "CH4") {
        return(ppt/1000*M_CH4()/Vd)
    } else {
        print("WARNING: You selected a gas which is not available for the conversion here.")
    }
    
}