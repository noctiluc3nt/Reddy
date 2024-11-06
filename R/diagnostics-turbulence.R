### turbulence quantities ###

#' Turbulent Kinetic Energy TKE
#'
#'@description Calculates turbulent kinetic energy (TKE) from \code{u_sd}, \code{v_sd} and \code{w_sd}
#'@param u_sd standard deviation of u-wind [m/s]
#'@param v_sd standard deviation of v-wind [m/s]
#'@param w_sd standard deviation of w-wind [m/s]
#'
#'@return turbulent kinetic energy TKE [m^2/s^2]
#'@export
#'
#'@examples
#'calc_tke(1,1,1)
#'
calc_tke = function(u_sd,v_sd,w_sd) {
	return(1/2*(u_sd^2+v_sd^2+w_sd^2)) 
}


#' Turbulent Kinetic Energy Velocity Scale
#'
#'@description Calculates the velocity scale of turbulent kinetic energy (TKE): \code{Vtke = sqrt(TKE)}
#'@param u_sd standard deviation of u-wind [m/s]
#'@param v_sd standard deviation of v-wind [m/s]
#'@param w_sd standard deviation of w-wind [m/s]
#'
#'@return turbulent kinetic energy velocity scale [m/s]
#'@export
#'
#'@examples
#'calc_vtke(1,1,1)
#'
calc_vtke = function(u_sd,v_sd,w_sd) {
	tke=calc_tke(u_sd,v_sd,w_sd)
	return(sqrt(tke)) 
}


#' Friction Velocity
#'
#'@description Calculates friction velocity from the covariances cov(u,w) and cov(v,w)
#'@param cov_uw covariance cov(u,w) [m^2/s^2]
#'@param cov_vw covariance cov(v,w) [m^2/s^2]
#'
#'@return friction velocity [m/s]
#'@export
#'
#'@examples
#'calc_ustar(1,1)
#'
calc_ustar = function(cov_uw,cov_vw) {
	return((cov_uw^2+cov_vw^2)^(1/4)) 
}

#' Obukhov length
#'
#'@description Calculates Obukhov length from friction velocity, mean temperature and cov(T,w)
#'@param ustar friction velocity (e.g., from \code{calc_ustar}) [m/s]
#'@param T_mean mean temperature [K]
#'@param cov_wT covariance cov(w,T) [m/s K]
#'
#'@return Obukhov length [m]
#'@export
#'
#'@examples
#'calc_L(0.2,273,0.1) #unstable
#'calc_L(0.2,273,-0.1) #stable
#'
calc_L = function(ustar,T_mean,cov_wT) {
	return(-abs(ustar^3)*T_mean/(karman()*g()*cov_wT))
}


#' Stability Parameter
#'
#'@description Calculates dimensionless stability parameter from Obukhov length and measurement height, i.e. \code{zeta = z/L}
#'@param z measurement height [m]
#'@param L Obukhov length [m] (e.g., from \code{calc_L})
#'
#'@return stability parameter [-]
#'@export
#'
#'@examples
#'calc_zeta(2,-1) #unstable
#'calc_zeta(2,1) #stable
#'
calc_zeta = function(z,L) {
	return(z/L)
}


#' Horizontal Turbulence Intensity TI
#'
#'@description Calculates horizontal turbulence intensity \code{TI = sqrt(u_sd^2+v_sd^2)/ws_mean}
#'@param u_sd standard deviation of streamwise wind (u-wind)
#'@param v_sd standard deviation of crosswise wind (v-wind)
#'@param ws_mean horizontal wind speed
#'
#'@return horizontal turbulence intensity [-]
#'@export
#'
#'@examples
#'calc_ti(1,1,3)
#'
calc_ti = function(u_sd,v_sd,ws_mean) {
	return(sqrt(u_sd^2+v_sd^2)/ws_mean)
}

#' Vertical Turbulence Intensity Iw
#'
#'@description Calculates vertical turbulence intensity \code{Iw = w_sd/ws_mean}
#'@param w_sd standard deviation of vertical wind (w-wind)
#'@param ws_mean horizontal wind speed
#'
#'@return vertical turbulence intensity [-]
#'@export
#'
#'@examples
#'calc_iw(1,3) #unstable
#'
calc_iw = function(w_sd,ws_mean) {
	return(w_sd/ws_mean)
}

#' Velocity Aspect Ratio (VAR)
#'
#'@description Calculates the velocity aspect ratio: \code{VAR = sqrt(2)*w_sd/sqrt(u_sd^2+v_sd^2)}
#'@param u_sd standard deviation of streamwise wind (u-wind)
#'@param v_sd standard deviation of crosswise wind (v-wind)
#'@param w_sd standard deviation of vertical wind (w-wind)
#'
#'@return velocity aspect ratio [-]
#'@export
#'
#'@examples
#'calc_var(1,1,1) #"isotropic"
#'calc_var(1,1,2) #not isotropic
#'
calc_var = function(u_sd,v_sd,w_sd) {
	return(sqrt(2)*w_sd/sqrt(u_sd^2+v_sd^2))
}

#' Directional Shear
#'
#'@description Calculates a measure for directional shear alpha_uw = arctan(cov(v,w)/cov(u,w))
#'@param cov_uw covariance cov(u,w)
#'@param cov_vw covariance cov(v,w)
#'
#'@return angle that describes the impact of directional shear [deg]
#'@export
#'
#'@examples
#'calc_dshear(-0.5,0) #no shear
#'calc_dshear(-0.5,-0.1)
#'
calc_dshear = function(cov_uw,cov_vw) {
	return((atan2(cov_vw,cov_uw)*180/pi+360)%%360)
}

#' Decoupling metric (Omega)
#'
#'@description Calculates the decoupling metric (Omega) from Peltola et al., 2021 (without vegetation)
#'@param w_sd standard deviation of vertical velocity [m/s]
#'@param N Brunt-Vaisala frequency [1/s]
#'@param z measurement height [m]
#'
#'@return decoupling metric (Omega) [-]
#'@export
#'
#'@examples
#'calc_decoupling_metric(1,1)
#'
calc_decoupling_metric = function(w_sd,N,z=2) {
	LB=w_sd/N #buoyancy length scale
	return(LB/(sqrt(2)*z)) #Peltola et al, 2021: eq 6
}

### intermittency indicator ###

#' Flux intermittency
#'
#'@description Calculates flux intermittency FI = flux_sd/flux (flux_sd: sd of subsampled fluxes) following Mahrt, 1998 (similar to stationarity flag \code{flag_stationarity})
#'@param ts1 timeseries 1 
#'@param ts2 timeseries 2 (optional), if the flux should be calculated based on \code{ts1*ts2} (default \code{ts2=NULL}, i.e. \code{ts2} is not used)
#'@param nsub number of elements used for subsampling, default \code{nsub=6000}, which corrosponds to 5 minutes of measurements from 20 Hz sampled half-hour (containing 30*60*20 = 36000 measurements)
#'
#'@return flux intermittency [-]
#'@export
#'
#'@examples
#'set.seed(5)
#'ts1=rnorm(30)
#'ts2=rnorm(30)
#'calc_flux_intermittency(ts1,ts2,nsub=6) #as product
#'calc_flux_intermittency(ts1*ts2,nsub=6) #the same from one variable
#'
calc_flux_intermittency = function(ts1,ts2=NULL,nsub=6000) {
	n=length(ts)
    nint=n%/%nsub
    if (nint<=1) {
        warning("nsub is chosen to large.")
    }
    cov_subs=array(NA,dim=nint)
    if (!is.null(ts2)) {
		cov_complete=cov(ts1,ts2,use="pairwise.complete.obs")
	} else {
		cov_complete=mean(ts1,na.rm=TRUE)
	}
    for (i in 1:nint) {
        isub=((i-1)*nsub+1):(i*nsub)
        if (!is.null(ts2)) cov_sub=cov(ts1[isub],ts1[isub],use="pairwise.complete.obs")
		if (is.null(ts2)) cov_sub=mean(ts1[isub],na.rm=TRUE)
        cov_subs[i]=cov_sub
    }
    return(sd(cov_subs)/cov_complete)
}


### hydrological measures ###
#' Bowen ratio BR
#'
#'@description Calculates the Bowen ratio as ratio of sensible and latent heat flux, i.e., BR := SH/LH
#'@param sh sensible heat flux [W/m^2]
#'@param lh latent heat flux [W/m^2]
#'
#'@return Bowen ratio [-]
#'@export
#'
#'@examples
#'calc_br(50,20)
#'
calc_br = function(sh,lh) {
	return(sh/abs(lh))
}

#' Evaporative fraction
#'
#'@description Calculates the evaporative fraction EF := LH/(SH+LH)
#'@param sh sensible heat flux [W/m^2]
#'@param lh latent heat flux [W/m^2]
#'
#'@return evaporative fraction [-]
#'@export
#'
#'@examples
#'calc_ef(50,20)
#'
calc_ef = function(sh,lh) {
	return(lh/(sh+lh))
}

#' Evapotranspiration
#'
#'@description Converts latent heat flux to evaporation
#'@param lh latent heat flux [W/m^2]
#'@param temp temperature [K] (optional), if provided, the latent heat of vaporization is calculated temperature-dependent
#'
#'@return evapotranspiration [kg/(s*m^2)]
#'@export
#'
#'@examples
#'lh2et(20)
#'lh2et(20,273)
#'
lh2et = function(lh,temp=NULL) {
	lv=Lv(temp)
	return(lh/lv)
}


### surface roughness and related concepts ###

#' Calculates surface roughness length z0 from friction velocity using the simple estimate from Charnock, 1955
#'
#'@description Calculates surface roughness z0 from friction velocity using the simple estimate from Charnock, 1955: z0 = alpha*ustar^2/g with alpha=0.016 and g=9.81 m/s^2
#'@param ustar friction velocity [m/s]
#'
#'@return surface roughness length [m]
#'@export
#'
#'@examples
#'ustar2z0(0.2)
#'
ustar2z0 = function(ustar) {
	return(alpha()*ustar^2/g())
}

