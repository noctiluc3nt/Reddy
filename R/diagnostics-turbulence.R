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
#'@examples-
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
#'@param cov_vw covariance cov(v,w) [m^2/s^2] (optional)
#'
#'@return friction velocity [m/s]
#'@export
#'
#'@examples
#'calc_ustar(-0.3,0.02)
#'
calc_ustar = function(cov_uw,cov_vw=0) {
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

#' Ozmidov scale (L_OZ)
#'
#'@description Calculates the Ozmidov length scale L_OZ = sqrt(epsilon/N^3), with epsilon: TKE dissipation rate, and N: Brunt-Vaisala frequency
#'@param epsilon dissipation rate of TKE or w [m*s]
#'@param N Brunt-Vaisala frequency [1/s]
#'
#'@return Ozmidov length scale [m]
#'@export
#'
#'@examples
#'calc_ozmidov_scale(-5/3,1*10^-4)
#'
calc_ozmidov_scale = function(epsilon,N) {
	return(sqrt(epsilon/N^3))
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
	n=length(ts1)
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



### Ekman layer

#' Coriolis parameter
#'
#'@description Calculates Coriolis parameter from latitude
#'@param phi latitude [deg]
#'
#'@return Coriolis parameter [1/s]
#'@export
#'
#'@examples
#'calc_coriolis(45)
#'
calc_coriolis = function(phi) {
	Omega=1/86400 #earth rotation frequency
	return(2*Omega*sin(phi*pi/180))
}

#' Ekman layer thickness
#'
#'@description Calculates Ekman layer thickness from eddy diffusivity and Coriolis parameter sqrt(2*Km/abs(f))
#'@param Km eddy diffusivity [m^2/s]
#'@param f Coriolis parameter [1/s] (e.g. from \code{calc_coriolis})
#'
#'@return Ekman layer thickness [m] (derived from boundary layer equations)
#'@export
#'
#'@examples
#'calc_ekman_layer_depth(0.1,10^(-4))
#'
calc_ekman_layer_depth = function(Km,f) {
	return(sqrt(2*Km/abs(f)))
}


### Brunt-Vaisala frequency + (bulk/flux) Richardson number

#' Brunt-Vaisala frequency squared
#'
#'@description calculates Brunt-Vaisala frequency squared (N^2)
#'@param T1 temperature at the lower level [K]
#'@param T2 temperature at the upper level [K]
#'@param dz height difference of the two measurements [m]
#'
#'@return N2 [1/s^2]
#'@export
#'
calc_N2 = function(T1,T2,dz) {
	T0=(T1+T2)/2
	dT_dz=(T2-T1)/dz
	return(T0/g()*dT_dz)
}

#' Calculates bulk Richardson number Ri
#'
#'@description calculates Richardson number Ri
#'@param U1 wind speed at the lower level [m/s]
#'@param U2 wind speed at the upper level [m/s]
#'@param T1 temperature at the lower level [K]
#'@param T2 temperature at the upper level [K]
#'@param dz height difference of the two measurements [m]
#'
#'@return Ri [-]
#'@export
#'
calc_ri = function(T1,T2,U1,U2,dz) {
	T0=(T1+T2)/2
	dT_dz=(T2-T1)/dz
	dUbar_dz=(U2-U1)/dz
	ri=g()/T0*dT_dz/(dUbar_dz^2)
	return(ri)
}

#' Calculates flux Richardson number Ri_f
#'
#'@description calculates flux Richardson number Ri_f = g/T_mean*cov(w,T)/(cov(u,w)*du/dz)
#'@param cov_wT covariance cov(w,T) [K m/s]
#'@param cov_uw covariance cov(u,w) [m^2/s^2]
#'@param U1 wind speed at the lower level [m/s]
#'@param U2 wind speed at the upper level [m/s]
#'@param dz height difference of the two measurements [m]
#'@param T_mean mean temperature [K] (optional, used instead of T0=273.15)
#'
#'@return Ri_f [-]
#'@export
#'
calc_rif = function(cov_wT,cov_uw,U1,U2,dz,T_mean=NULL) {
	T0=ifelse(is.null(T_mean),273.15,T_mean)
	dUbar_dz=(U2-U1)/dz
	rif=g()/T0*cov_wT/(cov_uw*dUbar_dz)
	return(rif)
}


### eddy viscosity and conductivity, Prandtl number

#' Calculates eddy viscosity K_m = -cov(u,w)/(du/dz)
#'
#'@description Calculates eddy viscosity K_m
#'@param cov_uw covariance cov(u,w) [m^2/s^2]
#'@param du_dz vertical wind speed gradient [1/s]
#'
#'@return eddy viscosity K_m [m^2/s]
#'@export
#'
#'@examples
#'calc_Km(-0.2,2)
#'
calc_Km = function(cov_uw,du_dz) {
	return(-cov_uw/du_dz)
}

#' Calculates eddy conductivity K_h = -cov(w,T)/(dT/dz)
#'
#'@description Calculates eddy conductivity K_h
#'@param cov_wT covariance cov(w,T) [K m/s]
#'@param dT_dz vertical temperature gradient [K/m]
#'
#'@return eddy conductivity K_h [m^2/s]
#'@export
#'
#'@examples
#'calc_Kh(0.2,-1)
#'
calc_Kh = function(cov_wT,dT_dz) {
	return(-cov_wT/dT_dz)
}

#' Calculates turbulent Prandtl number Pr = K_m/K_h
#'
#'@description Calculates turbulent Prandtl number Pr
#'@param K_m eddy viscosity [m^2/s]
#'@param K_h eddy conductivity [m^2/s]
#'
#'@return Prandtl number [-]
#'@export
#'
#'@examples
#'calc_Pr(0.4,0.6)
#'
calc_Pr = function(K_m,K_h) {
	return(K_m/K_h)
}