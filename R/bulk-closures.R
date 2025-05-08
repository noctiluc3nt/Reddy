### scaling functions used for flux-variance relations ###

#' Scaling function for horizontal windspeed Phi_u
#'
#'@description scaling function Phi_u 
#'@param zeta stability parameter [-]
#'@param method defining from which paper the scaling function is used, default \code{method="PD1984"} for Panofsky and Dutton, 1984
#'
#'@return Phi_u
#'@export
#'
scale_phiu = function(zeta,method="PD1984") {
    if (method=="PD1984") {
        if (zeta<=0) { #unstable
            return(2.55*(1-3*zeta)^(1/3))
        } else { #stable
            return(2.55*(1+3*zeta)^(1/3))
        }
    }
}

#' Scaling function for vertical windspeed Phi_w
#'
#'@description scaling function Phi_w
#'@param zeta stability parameter [-]
#'@param method defining from which paper the scaling function is used, default \code{method="PD1984"} for Panofsky and Dutton, 1984
#'
#'@return Phi_w
#'@export
#'
scale_phiw = function(zeta,method="PD1984") {
    if (method=="PD1984") {
	    if (zeta<=0) { #unstable
            return(1.25*(1-3*zeta)^(1/3))
        } else { #stable
          return(1.25*(1+3*zeta)^(1/3))
        }
    }
}

#' Scaling function for temperature Phi_T
#'
#'@description scaling function Phi_T
#'@param zeta stability parameter [-]
#'@param method defining from which paper the scaling function should be used, default \code{method="K1994"} for Katul, 1994, other option \code{method="SC2018"} for Stiperski and Calaf, 2018
#'
#'@return Phi_T
#'@export
#'
scale_phiT = function(zeta,method="K1994") {
    if (method=="K1994") {
        if (zeta<=0) { #unstable
            return(0.95*(-zeta)^(-1/3))
        } else { #stable
            return(0.95*(zeta)^(-1/3))
        }
    } else if (method=="SC2018") {
        if (zeta<=0) { #unstable
            return(0.99*(0.067-zeta)^(-1/3))
        } else { #stable
            return(1.76+0.15*(zeta)^(-1))
        }
    }	
}



### Scaling functions for flux-profile relations ###

#' Scaling function for momentum Phi_m
#'
#'@description scaling function Phi_m
#'@param zeta stability parameter [-]
#'@param method defining from which paper the scaling function should be used, default \code{method="ecmwf"} for the ones used in ECMWF-IFS, other option \code{method="BD"} for the linear Businger-Dyer relations
#'
#'@return Phi_m
#'@export
#'
scale_phim = function(zeta,method="ecmwf") {
    if (method=="ecmwf") {
        if (zeta<=0) { #unstable
            return((1-16*zeta)^(-1/4))
        } else { #stable
            return(1+5*zeta)
        }
    }
    if (method=="BD") {
        if (zeta<=0) { #unstable
            return(2-zeta)
        } else { #stable
            return(1+5*zeta)
        }
    }	
}

#' Scaling function for heat Phi_h
#'
#'@description scaling function Phi_h
#'@param zeta stability parameter [-]
#'@param method defining from which paper the scaling function should be used, default \code{method="ecmwf"} for for the ones used in ECMWF-IFS
#'
#'@return Phi_h
#'@export
#'
scale_phih = function(zeta,method="BD") {
    if (method=="ecmwf") {
        if (zeta<=0) { #unstable
            return((1-16*zeta)^(-1/2))
        } else { #stable
            return((1+4*zeta)^2)
        }
    }
}



### calculate scaling parameters for flux-profile relations ###

#' Calculates Phi_m
#'
#'@description calculates scaling function Phi_m (for momentum)
#'@param U1 wind speed at the lower level [m/s]
#'@param U2 wind speed at the upper level [m/s]
#'@param ustar friction velocity [m/s]
#'@param zm measurement/scaling height [m]
#'@param dz height difference of the two measurements [m]
#'
#'@return Phi_m
#'@export
#'
calc_phim = function(U1,U2,ustar,zm,dz) {
	dUbar_dz=(U2-U1)/dz
	phi=karman()*zm*abs(dUbar_dz)/ustar
	return(phi)
}

#' Calculates Phi_h
#'
#'@description calculate scaling function Phi_h (for heat)
#'@param T1 temperature at the lower level [K]
#'@param T2 temperature at the upper level [K]
#'@param cov_wT covariance cov(w,T) [K m/s]
#'@param ustar friction velocity [m/s]
#'@param zm measurement/scaling height [m]
#'@param dz height difference of the two measurements [m]
#'
#'@return Phi_h
#'@export
#'
calc_phih = function(T1,T2,cov_wT,ustar,zm,dz) {
    tstar=calc_xstar(cov_wT,ustar)
	dTbar_dz=(T2-T1)/dz
	phi=karman()*zm*abs(dTbar_dz)/tstar
	return(phi)
}


#' Calculates xstar (denominator for general flux-variance relation)
#'
#'@description calculates xstar = x/ustar (for general flux-variance relation)
#'@param x variable that should be scaled
#'@param ustar friction velocity [m/s]
#'
#'@return xstar = x/ustar
#'@export
#'
calc_xstar = function(x,ustar) {
	return(x/ustar)
}


#' Calculates Phi_x (general flux-variance relation)
#'
#'@description calculates Phi_x = sigma_x/xstar (for general flux-variance relation)
#'@param sigma_x standard deviation of x
#'@param x variable that should be scaled, e.g. vertical flux of x with x = T or x = q
#'@param ustar friction velocity [m/s]
#'
#'@return Phi_x = sigma_x/xstar
#'@export
#'
calc_phix = function(sigma_x,x,ustar) {
    xstar=calc_xstar(x,ustar)
	return(sigma_x/xstar)
}


#' Wind profile from Monin-Obukhov similarity theory
#'
#'@description Calculates vertical profile of horizontal wind speed following Monin-Obukhov similarity theory
#'@param zs scalar or vector, heights [m] at which the horizontal wind speed should be calculate 
#'@param ustar friction velocity [m/s]
#'@param z0 surface roughness length [m], default \code{z0=0} (note: it could be an option to calculate z0 from ustar with \code{ustar2z0()})
#'@param d displacement height [m], optional, default \code{d=0} (i.e. no displacement)
#'@param zeta stability parameter [-] to correct for stability effects, default \code{zeta=0} (i.e. no stability correction, resulting in classical logarithmic wind profile)
#'@param method "method" for calculating stability correction function (only relevant if zeta is non-zero), default \code{method="ecmwf"} for using Phi_m from ECMWF-IFS
#'
#'@return data frame containing the requested heights \code{zs} and the calculated wind speed [m/s] there
#'@export
#'
#'@examples
#'zs=seq(1,100)
#'ustar=0.2
#'u_neutral=calc_windprofile(zs,ustar)
#'u_unstable=calc_windprofile(zs,ustar,zeta=-0.2)
#'u_stable=calc_windprofile(zs,ustar,zeta=0.2)
#'
calc_windprofile = function(zs,ustar,z0=0,d=0,zeta=0,method="ecmwf") {
    Phi=scale_phim(zeta,method)
    if (z0==0) {
        uz=ustar/karman()*(log(zs-d)+Phi)
    } else {
        uz=ustar/karman()*log((zs-d)/z0+Phi)
    }
	return(data.frame("height"=zs,"windspeed"=uz))
}
