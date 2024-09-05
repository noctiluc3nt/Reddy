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
#'@param method defining from which paper the scaling function should be used, default \code{method="ecmwf"} for for the ones used in ECMWF-IFS, other option \code{method="BD"} for the lienar Businger-Dyer relations
#'
#'@return Phi_m
#'@export
#'
scale_phim = function(zeta,method="BD") {
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