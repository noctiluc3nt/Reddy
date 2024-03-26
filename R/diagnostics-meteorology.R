#' Saturation vapor pressure over water
#'
#'@description Calculates the saturation vapor pressure over water for given temperature and pressure
#'@param temp scalar or vector, temperature [°C]
#'@return E_s, saturation vapor pressure over water [hPa]
#'@export
#'
calc_satvaporpressure = function(temp) {
    a=0.61094
    b=17.625
    c=243.04
    return(a*exp(b*temp/(temp+c))/10)
}


#' Clear Sky Index (CSI)
#'
#'@description Calculates clear sky index
#'@param temp scalar or vector, temperature [K]
#'@param lw_in scalar or vector, longwave incoming radiation [W/m^2]
#'@param rh scalar or vector, relative humidity [percent]
#'@param e scalar or vector, vapor pressure [Pa]
#'@return CSI, clear sky index
#'@export
#'
calc_csi = function(temp,lw_in,rh=NULL,e=NULL) {
    if (is.null(rh) & is.null(e)) {
        print("ERROR: Either relative humidity rh or vapor pressure e have to be given.")
    }
    if (!is.null(rh)) { #calculate vapor pressure
        es = calc_satvaporpressure(temp-273.15)
        e = rh * es/100
    }
    sigma=5.67*10^(-8)
    epsilon_A = lw_in/(sigma*temp^4) #actual atmospheric emissivity
    epsilon = 0.23 + 0.47*(100*e/temp)^(1/8) #(theoretical) clear sky emissivity, Marty and Philiponna, 2002
    return(epsilon_A/epsilon)
}