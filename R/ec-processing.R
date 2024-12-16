#' Despiking
#'
#'@description Applies (up to) three despiking methods based on (1) pre-defined thresholds, (2) median deviation (MAD) test and (3) skewness and kurtosis
#'@param ts timeseries that shall be despiked
#'@param thresholds vector with two elements representing lower and upper bounds for despiking (pre-defined thresholds), \code{NA} means that the respective bound is not used
#'@param mad_factor factor for the MAD test, default \code{mad_factor = 10}
#'@param threshold_skewness threshold for skewness test, default \code{threshold_skewness = 2}
#'@param threshold_kurtosis threshold for kurtosis test, default \code{threshold_kurtosis = 8}
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
despiking = function(ts,thresholds=c(NA,NA),mad_factor=10,threshold_skewness=2,threshold_kurtosis=8) {
	#despiking based on predefined limits
    if (sum(is.na(thresholds))==0) {
        pass = (ts>thresholds[1] & ts<thresholds[2])
	    ts[!pass] = NA
    } else if (is.na(thresholds[1])) {
        pass = (ts<thresholds[2])
	    ts[!pass] = NA
    } else if (is.na(thresholds[2])) {
        pass = (ts>thresholds[1])
	    ts[!pass] = NA
    }
    #despiking based on median deviation test
    med=median(ts,na.rm=TRUE)
    mad=median(abs(ts-med),na.rm=TRUE)
    pass=(abs(ts-med) <= mad_factor*mad) #pass criterion
    ts[!pass] = NA
    #despiking based on skewness and kurtosis for entire block
    seriesLDT=pracma::detrend(ts,tt="linear") #linear detrending to eliminate trends (departures from stationarity) -> would influnece higher moments
    skewness=mean(seriesLDT^3,na.rm=T)/sd(seriesLDT,na.rm=T)^3
    kurtosis=mean(seriesLDT^4,na.rm=T)/sd(seriesLDT,na.rm=T)^4
    pass=(abs(skewness)<threshold_skewness & kurtosis<threshold_kurtosis)
    if (pass==FALSE) ts = array(NA,dim=length(ts))
    return(ts)
}


#' Double rotation
#'
#'@description Double rotation (i.e., sonic coordinate system will be aligned with streamlines)
#'@param u u-wind (levelled sonic)
#'@param v v-wind (levelled sonic)
#'@param w w-wind (levelled sonic)
#'
#'@return list containing the wind in a natural coordinate system (streamwise, crosswise, vertical) and the two rotation angles theta and phi
#'@export
#'
#'@examples
#'wind_rotated=rotate_double(4,3,1) #double rotation can be applied instantenously
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
    return(list("u"=u2,"v"=v2,"w"=w2,"theta"=(theta*180/pi+360)%%360,"phi"=(phi*180/pi+360)%%360))
}


#' Planar fit rotation
#'
#'@description Planar fit rotation (i.e., sonic coordinate system will be aligned with the mean streamlines resulting in vanishing of w_mean) 
#'@param u u-wind (levelled sonic)
#'@param v v-wind (levelled sonic)
#'@param w w-wind (levelled sonic)
#'@param bias a three-dimensional correction vector containing the offsets of u-, v-, w-wind
#'
#'@return list containing u, v, w after planar fit rotation as well as the rotation angles alpha, beta and gamma and the fitted offset c3
#'@export
#'
#'@examples
#'u=rnorm(1000)
#'v=rnorm(1000)
#'w=rnorm(1000)
#'wind_rotated=rotate_planar(u,v,w) #for planar fit a timeseries is required
#'
rotate_planar = function(u,v,w,bias=c(0,0,0)) {
    if (!identical(length(u),length(v),length(w))) { 
        stop("u, v, w have to be of same length.")
    }
    #linear regression
    fit=lm(w ~ u + v)
    c3=fit$coefficients[1]
    b1=fit$coefficients[2]
    b2=fit$coefficients[3]
    #angles
    alpha=atan(-b1)
    beta=atan(b2)
    sin_alpha=-b1/(sqrt(1+b1^2))
    cos_alpha=1/(sqrt(1+b1^2))
    sin_beta=b2/(sqrt(1+b2))
    cos_beta=1/(sqrt(1+b2^2))
    #rotation matrix P
    P=matrix(c(cos_alpha,sin_alpha*sin_beta,-sin_alpha*cos_beta,0,cos_beta,sin_beta,sin_alpha,-sin_beta*cos_alpha,cos_alpha*cos_beta),nrow=3,byrow=T)
    upf=P[1,1]*(u-bias[1]) + P[1,2]*(v-bias[2]) + P[1,3]*(w-bias[3])
    vpf=P[2,1]*(u-bias[1]) + P[2,2]*(v-bias[2]) + P[2,3]*(w-bias[3])
    wpf=P[3,1]*(u-bias[1]) + P[3,2]*(v-bias[2]) + P[3,3]*(w-bias[3])
    #gamma angle for rotation around z-axis
    gamma=atan2(mean(vpf,na.rm=T),mean(upf,na.rm=T))
    ur=upf*cos(gamma) + vpf*sin(gamma)
    vr=-upf*sin(gamma) + vpf*cos(gamma)
    return(list("u"=ur,"v"=vr,"w"=wpf,"alpha"=(alpha*180/pi+360)%%360,"beta"=(beta*180/pi+360)%%360,"gamma"=(gamma*180/pi+360)%%360,"c3"=c3))
}


#' Stationarity Flag
#'
#'@description Stationarity Flag according to Foken and Wichura, 1996 based on the assumption that the covariance of two variables (\code{var1} and \code{var2}, one usually representing vertical velocity) calculated for blocks (of length \code{nsub}) does not differ to much from the total covariance
#'@param var1 variable 1 
#'@param var2 variable 2 (same length as \code{var1}, usually either \code{var1} or \code{var2} represent vertical velocity)
#'@param nsub number of elements used for subsampling (\code{nsub < length(var1)}) 
#'@param thresholds_stationarity vector containing 2 elements to distinguish between flag=0 and flag=1, as well as flag=1 and flag=2, default: \code{c(0.3,1)}
#'
#'@return stationarity flags (0: in full agreement with the criterion ... 2: does not fulfill the criterion)
#'@export
#'
#'@examples
#'set.seed(5)
#'ts1=rnorm(30)
#'ts2=rnorm(30)
#'flag_stationarity(ts1,ts2,nsub=6)
#'
flag_stationarity = function(var1,var2,nsub=3000,thresholds_stationarity=c(0.3,1)) {
    if (length(var1) != length(var2)) {
        stop("var1 and var2 have to be of equal length.")
    }
    if (length(thresholds_stationarity)!=2) {
        stop("thresholds_stationarity has to be a vector of length 2.")
    }
    nint=length(var1)%/%nsub
    if (nint<=1) {
        stop("nsub is chosen to large.")
    }
    rnfs=array(NA,dim=nint)
    cov_complete=cov(var1,var2,use="pairwise.complete.obs")
    for (i in 1:nint) {
        isub=((i-1)*nsub+1):(i*nsub)
        cov_sub=cov(var1[isub],var2[isub],use="pairwise.complete.obs")
        rnfs[i]=abs((cov_complete-cov_sub)/cov_complete)
    }
    rnf=mean(rnfs,na.rm=T)
    flag=ifelse(rnf<thresholds_stationarity[1],0,ifelse(rnf<thresholds_stationarity[2],1,2))
    return(flag)
}


#' Vertical Velocity Flag
#'
#'@description Vertical velocity flag according to Mauder et al., 2013: After rotation the vertical velocity should vanish, this flag flags high remaining vertical velocities.
#'@param w vertical velocity
#'@param thresholds_w vector containing 2 elements to distinguish between flag=0 and flag=1, as well as flag=1 and flag=2, default: \code{c(0.1,0.15)}
#'
#'@return vertical velocity flags (0: in full agreement with the criterion ... 2: does not fulfill the criterion)
#'@export
#'
#'@examples
#'flag_w(0.01)
#'
flag_w = function(w,thresholds_w=c(0.1,0.15)) {
    if (length(thresholds_w)!=2) {
        stop("thresholds_w has to be a vector of length 2.")
    }
    w=abs(w)
    flag=ifelse(w<thresholds_w[1],0,ifelse(w<thresholds_w[2],1,2))
    return(flag)
}

#' Flow Distortion Flag and Wind Constancy Ratio
#'
#'@description Flow Distortion Flag according to Mauder et al., 2013: Wind coming from (pre-defined) directions blocked by the measurement device is flaged with 2 (for wind speeds greater than 0.1 assuming that during calm wind the wind direction is not well-defined). The wind constancy ratio is calculated to quantify the variability of horizontal wind direction according to Mahrt, 1999.
#'@param u u-wind (levelled sonic)
#'@param v v-wind (levelled sonic)
#'@param dir_blocked vector containing the lower and upper bound of the blocked wind sector in degrees (e.g., \code{dir_blocked = c(30,60)})
#'@param threshold_cr threshold for constancy ratio (default \code{threshold_cr = 0.9}, may be adapted to used data set)
#'
#'@return distortion flags (0: in full agreement with the criterion ... 2: does not fulfill the criterion)
#'@export
#'
#'@examples
#'flag_distortion(1,1,dir_blocked=c(30,60))
#'flag_distortion(1,1,dir_blocked=c(180,360))
#'
flag_distortion = function(u,v,dir_blocked,threshold_cr=0.9) {
    if (length(u) != length(v)) {
        stop("u and v have to be of equal length.")
    }
    #horizontal wind speed
    ws=sqrt(mean(u,na.rm=TRUE)^2+mean(v,na.rm=TRUE)^2)
    #constancy ratio cr
    cr=sqrt(sum(u^2,na.rm=TRUE)+sum(v^2,na.rm=TRUE))/(ws)
    #flow distortion flag considering cr
    if (!is.na(ws) & !is.na(cr)) {
        if (ws>0.1 & cr>threshold_cr) {
            wd=median(calc_windDirection(u,v))
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
#'@description Integral Turbulence Characteristics Flag: Tests the consistency with Monin-Obukhov similarity theory using the scaling functions from Panofsky and Dutton, 1984.
#'@param w_sd standard deviation of vertical velocity
#'@param ustar friction velocity
#'@param zeta stability parameter \code{zeta = z/L}
#'@param thresholds_most vector containing 2 elements to distinguish between flag=0 and flag=1, as well as flag=1 and flag=2, default: \code{c(0.3,0.8)}
#'
#'@return integral turbulence characteristics flags (0: in full agreement with the criterion ... 2: does not fulfill the criterion)
#'@export
#'
#'@examples
#'itc_flag=flag_most(0.2,0.4,-0.3)
#'
flag_most = function(w_sd,ustar,zeta,thresholds_most=c(0.3,0.8)) {
    if (length(thresholds_most)!=2) {
        stop("thresholds_most has to be a vector of length 2.")
    }
    parameterized=1.3*(1+2*abs(zeta))^(1/3) #w_sd/ustar parametrized according to scaling function based on zeta
    itc=abs((w_sd/ustar-parameterized)/parameterized)
    flag=ifelse(itc<thresholds_most[1],0,ifelse(itc<thresholds_most[2],1,2))
    return(flag)
}

#' Ts2T
#'
#'@description Converts sonic temperature Ts to temperature T
#'
#'@param Ts sonic temperature [K] (similar as virtual temperature)
#'@param q specific humidity [kg/kg]
#'
#'@return temperature [K]
#'@export
#'
Ts2T = function(Ts,q) {  
    return(Ts*(1+Rd()/Rv()*q))
}

#' SND and cross-wind correction of sensible heat flux
#'
#'@description SND and cross-wind correction of sensible heat flux: converts the buoyancy flux cov(w,Ts) (based on sonic temperature Ts) to sensible heat flux
#'@param Ts_mean sonic temperature [K] (averaged)
#'@param u_mean u-wind [m/s] (averaged)
#'@param v_mean v-wind [m/s] (averaged)
#'@param cov_uw cov(u,w) [m^2/s^2]
#'@param cov_vw cov(v,w) [m^2/s^2]
#'@param cov_wTs cov(Ts,w) [K*m/s] (buoyancy flux)
#'@param cov_qw cov(q,w) [kg/kg*m/s] (optional)
#'@param A constant used in cross-wind correction, default \code{A = 7/8} for CSAT3
#'@param B constant used in cross-wind correction, default \code{B = 7/8} for CSAT3
#'@param sos speed of sound [m/s], default \code{sos = csound()} corresponding to 343 m/s
#'
#'@return SND correction of sensible heat flux
#'@export
#'
SNDcorrection = function(Ts_mean,u_mean,v_mean,cov_uw,cov_vw,cov_wTs,cov_qw=NULL,A=7/8,B=7/8,sos=csound()) {
    if (!is.null(cov_qw)) { #considering q
        #second term: SND correction, third term: cross-wind correction
        return(cov_wTs - 0.51*cov_qw + 2*Ts_mean/csound()^2*(A*u_mean*cov_uw + B*v_mean*cov_vw))
    }
    #without q: only cross-wind correction
    return(cov_wTs + 2*Ts_mean/csound()^2*(A*u_mean*cov_uw + B*v_mean*cov_vw))
}

#' WPL correction
#'
#'@description WPL correction: density correction for trace gas fluxes (i.e., converts volume- to mass-related quantity)
#'@param Ts_mean temperature [K] (sonic temperature or corrected temperature)
#'@param q_mean specific humidity [kg/kg] (if measured, default \code{NULL})
#'@param cov_wTs covariance cov(w,Ts) [m/s*K]
#'@param rhow_mean measured water vapor density [kg/m^3]
#'@param cov_wrhow covariance cov (w,rhow) [m/s*kg/m^3]
#'@param rhoc_mean measured trace gas density [kg/m^3] (only if WPL-correction should be applied to another flux, e.g. CO2 flux, default \code{NULL})
#'@param cov_wrhoc covariance cov (w,rhoc) [m/s*kg/m^3] (only if WPL-correction should be applied to another flux, e.g. CO2 flux, default \code{NULL})
#'
#'@return WPL correction of respective flux
#'@export
#'
WPLcorrection = function(Ts_mean,q_mean,cov_wTs,rhow_mean,cov_wrhow,rhoc_mean=NULL,cov_wrhoc=NULL) {
    if (is.null(rho_c)) { #water vapor flux
        return((1+1.61*q_mean)*(cov_wrhow+rhow_mean/Ts_mean*cov_wTs)) #with M_L/M_w = 1.61
    } else { #other trace gas flux
        return(cov_wrhoc+1.61*rhoc_mean/rhow_mean*cov_wrhow+(1+1.61*q_mean)*rhoc_mean/Ts_mean*cov_wTs)
    }
}


#' Unit conversion of "parts-per" to density (for closed-path gas analyzer)
#'
#'@description Unit conversion of "parts-per" to density (for closed-path gas analyzer)
#'@param ppt measurement in parts per thousand [ppt]
#'@param T_mean temperature [K]
#'@param pres pressure [Pa]
#'@param e water vapor pressure [Pa]
#'@param gas which gas? can be either \code{H2O}, \code{CO2}, \code{CH4} (if \code{CO2}/\code{CH4} is selected, make sure that it's still in ppt and not ppm as usual)
#'
#'@return density of the gas [kg/m^3]
#'@export
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
        warning("You selected a gas which is not available for the conversion here.")
    } 
}


### speed of sound to sonic temperature ###
#' Converts speed of sound (sos) to sonic temperature
#'
#'@description Converts speed of sound (sos) to sonic temperature
#'@param sos speed of sound [m/s]
#'
#'@return sonic temperature (virtual temperature) [K]
#'@export
#'
sos2Ts = function(sos) {
	return(sos^2/(cpcv()*Rd()))
}


### fluxes (unit conversion) ###
#' Converts cov(w,T) to sensible heat flux SH
#'
#'@description Converts cov(T,w) to sensible heat flux SH
#'@param cov_wT covariance cov(w,T) [K m/s]
#'@param rho density of air [kg/m^3] (optional)
#'
#'@return sensible heat flux [W/m^2]
#'@export
#'
cov2sh = function(cov_wT,rho=NULL) {
	if (is.null(rho)) {
		rho = rhoAir()
	}
	return(rho*cp()*cov_wT)
}

#' Converts cov(w,q) to latent heat flux LH
#'
#'@description Converts cov(w,q) to latent heat flux LH
#'@param cov_wq covariance cov(w,q) [m/s]
#'@param rho density of air [kg/m^3] (optional)
#'
#'@return latent heat flux [W/m^2]
#'@export
#'
cov2lh = function(cov_wq,rho=NULL) {
	if (is.null(rho)) {
		rho = rhoAir()
	}
	return(rho*Lv()*cov_wq)
}

#' Converts cov(co2,w) to CO2 flux
#'
#'@description Converts cov(co2,w) to CO2 flux
#'@param cov_co2w covariance cov(co2,w) [m/s]
#'@param rho density of air [kg/m^3] (optional)
#'
#'@return CO2 flux [kg/(m^2*s)]
#'@export
#'
cov2cf = function(cov_co2w,rho=NULL) {
	if (is.null(rho)) {
		rho = rhoAir()
	}
	return(rho*cov_co2w)
}


#' Calculates covariance of two timeseries using pair-wise complete observations
#'
#'@description Calculates cov(x,y)
#'@param x timeseries 1
#'@param y timeseries 2
#'
#'@return cov(x,y)
#'@export
#'
#'@examples
#'set.seed(5)
#'x=rnorm(100)
#'y=rnorm(100)
#'y[1:10]=NA
#'cov_xy=calc_cov(x,y)
#'
calc_cov = function(x,y) {
	return(cov(x,y,use="pairwise.complete.obs"))
}


#' Response-time correction factor (spectral correction)
#'
#'@description Calculates the response-time correction factor from cospectrum, e.g. Peltola et al., 2021
#'@param cospectrum cospectrum
#'@param freq frequency [Hz], corresponding to the cospectrum, i.e. same length
#'@param tau response time of the instrument [s] (has to be determined first by comparison with another instrument, that samples faster)
#'
#'@return response-time correction factor (which then can be used to correct the covariance and fluxes by multiplication)
#'@export
#'
RTcorrection = function(cospectrum,freq,tau=1) {
    tf=(1+(2*pi*tau*freq)^2) #transfer function that accounts for high-frequency loss due to limited sampling frequency
    cospectrum_cor=cospectrum*sqrt(tf) #note: their is a discussion whether sqrt(tf) or tf should be applied in the correction
    rt_factor=sum(cospectrum_cor)/sum(cospectrum)
	return(rt_factor)
}
