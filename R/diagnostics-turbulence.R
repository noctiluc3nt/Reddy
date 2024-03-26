### turbulence quantities ###

#' Turbulent Kinetic Energy
#'
#'@description Calculates turbulent kinetic energy (TKE) from u_sd, v_sd and w_sd
#'@param u_sd standard deviation of u-wind [m/s]
#'@param v_sd standard deviation of v-wind [m/s]
#'@param w_sd standard deviation of w-wind [m/s]
#'
#'@return turbulent kinetic energy TKE [m^2/s^2]
#'@export
#'
#'@examples
#'
calc_tke = function(u_sd,v_sd,w_sd) {
	return(1/2*(u_sd^2+v_sd^2+w_sd^2)) 
}


#' Turbulent Kinetic Energy Velocity Scale
#'
#'@description Calculates the velocity scale of turbulent kinetic energy (TKE): Vtke := sqrt(TKE)
#'@param u_sd standard deviation of u-wind [m/s]
#'@param v_sd standard deviation of v-wind [m/s]
#'@param w_sd standard deviation of w-wind [m/s]
#'
#'@return turbulent kinetic energy TKE [m^2/s^2]
#'@export
#'
#'@examples
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
#'
calc_ustar = function(cov_uw,cov_vw) {
	return((cov_uw^2+cov_vw^2)^(1/4)) 
}

#' Obukhov length
#'
#'@description Calculates Obukhov length from friction velocity, mean temperature and cov(T,w)
#'@param ustar friction velocity (e.g., from calc_frictionVelocity) [m/s]
#'@param T_mean mean temperature [K]
#'@param cov_wT covariance cov(w,T) [m/s K]
#'
#'@return Obukhov length [m]
#'@export
#'
#'@examples
#'
calc_L = function(ustar,T_mean,cov_wT) {
	return(-abs(ustar^3)*T_mean/(kap()*g()*cov_wT))
}


#' Stability Parameter
#'
#'@description Calculates dimensionless stability parameter from Obukhov length and measurement height, i.e. zeta = z/L
#'@param z measurement height [m]
#'@param L Obukhov length [m], e.g. from calc_L
#'
#'@return stability parameter [-]
#'@export
#'
#'@examples
#'
calc_zeta = function(z,L) {
	return(z/L)
}


#' Horizontal Turbulence Intensity TI
#'
#'@description Calculates horizontal turbulence intensity TI := sqrt(u_sd^2+v_sd^2)/ws_mean
#'@param u_sd standard deviation of streamwise wind (u-wind)
#'@param v_sd standard deviation of crosswise wind (v-wind)
#'@param ws_mean horizontal wind speed
#'
#'@return horizontal turbulence intensity [-]
#'@export
#'
#'@examples
#'
calc_ti = function(u_sd,v_sd,ws_mean) {
	return(sqrt(u_sd^2+v_sd^2)/ws_mean)
}

#' Vertical Turbulence Intensity Iw
#'
#'@description Calculates vertical turbulence intensity Iw := w_sd/ws_mean
#'@param w_sd standard deviation of vertical wind (w-wind)
#'@param ws_mean horizontal wind speed
#'
#'@return vertical turbulence intensity [-]
#'@export
#'
#'@examples
#'
calc_iw = function(w_sd,ws_mean) {
	return(w_sd/ws_mean)
}

#' Velocity Aspect Ratio (VAR)
#'
#'@description Calculates the velocity aspect ratio: VAR := sqrt(2)*w_sd/sqrt(u_sd^2+v_sd^2)
#'@param u_sd standard deviation of streamwise wind (u-wind)
#'@param v_sd standard deviation of crosswise wind (v-wind)
#'@param w_sd standard deviation of vertical wind (w-wind)
#'
#'@return velocity aspect ratio [-]
#'@export
#'
#'@examples
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
	return(atan2(cov_vw,cov_uw)*180/pi)
}


### fluxes ###
#' Converts cov(w,T) to sensible heat flux SH
#'
#'@description Converts cov(T,w) to sensible heat flux SH
#'@param cov_wT covariance cov(w,T) [K m/s]
#'@param rho density of air [kg/m^3] (optional)
#'
#'@return sensible heat flux [W/m^2]
#'@export
#'
cov2sh = function(covar_wT,rho=NULL) {
	if (is.null(rho)) {
		rho = rhoAir()
	}
	return(rho*cp()*covar_wT)
}

#' Converts cov(w,q) to latent heat flux LH
#'
#'@description Converts cov(q,w) to latent heat flux LH
#'@param cov_wT covariance cov(w,q) [m/s]
#'@param rho density of air [kg/m^3] (optional)
#'
#'@return latent heat flux [W/m^2]
#'@export
#'
cov2lh = function(covar_wq,rho=NULL) {
	if (is.null(rho)) {
		rho = rhoAir()
	}
	return(rho*covar_wq)
}

### wind basics ###
#' Wind Direction
#'
#'@description Calculates (horizontal) wind direction
#'@param u u-wind [m/s]
#'@param v v-wind [m/s]
#'
#'@return wind direction [deg]
#'@export
#'
#'@examples
#'
calc_windDirection = function(u,v) {
	return(atan2(-u,-v)*180/pi)
}

#' Horizontal Wind Speed
#'
#'@description Calculates horizontal wind speed
#'@param u u-wind [m/s]
#'@param v v-wind [m/s]
#'
#'@return wind speed [m/s]
#'@export
#'
#'@examples
#'
calc_windSpeed2D = function(u,v) {
	return(sqrt(u^2+v^2))
}

#' Wind Speed (3D)
#'
#'@description Calculates wind speed (3D)
#'@param u u-wind [m/s]
#'@param v v-wind [m/s]
#'@param w w-wind [m/s]
#'
#'@return wind speed (3D) [m/s]
#'@export
#'
#'@examples
#'
calc_windSpeed3D = function(u,v,w) {
	return(sqrt(u^2+v^2+w^2))
}

#' Gust Factor
#'
#'@description Calculates gust factor G := ws_max/ws_mean
#'@param ws_max wind speed [m/s]
#'@param ws_mean wind speed maximum [m/s]
#'
#'@return gust factor [-]
#'@export
#'
#'@examples
#'
calc_gustfactor = function(ws_max,ws_mean) {
	return(ws_max/ws_mean)
}