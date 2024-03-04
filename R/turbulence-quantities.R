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


#' Friction Velocity
#'
#'@description Calculates friction velocity from the covariances cov(u,w) and cov(v,w)
#'@param covar_uw covariance cov(u,w) [m^2/s^2]
#'@param covar_vw covariance cov(v,w) [m^2/s^2]
#'
#'@return friction velocity [m/s]
#'@export
#'
#'@examples
#'
calc_frictionVelocity = function(covar_uw,covar_vw) {
	return((covar_uw^2+covar_vw^2)^(1/4)) 
}

#' Obukhov length
#'
#'@description Calculates Obukhov length from friction velocity, mean temperature and cov(T,w)
#'@param ustar friction velocity (e.g., from calc_frictionVelocity) [m/s]
#'@param T_mean mean temperature [K]
#'@param covar_wT covariance cov(w,T) [m/s K]
#'
#'@return Obukhov length [m]
#'@export
#'
#'@examples
#'
calc_L = function(ustar,T_mean,covar_wT) {
	return(-abs(ustar^3)*T_mean/(kap()*g()*covar_wT))
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



### fluxes ###
#calc_Hflux = function(covar_wT,w_mean,T_mean,rho=NULL) {
#	if (is.null(rho)) {
#		rho = rhoAir()
#	}
	#return(rho*cp()*covar_wT+w_mean*T_mean)
#	return(rho*cp()*covar_wT)
#}

#calc_Eflux = function(covar_wq,w_mean,T_mean,rho=NULL) {
#	if (is.null(rho)) {
#		rho = rhoAir()
#	}
#	return(rho*covar_wq)
#}

### wind basic ###

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