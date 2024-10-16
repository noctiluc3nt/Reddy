#' Saturation vapor pressure over water
#'
#'@description Calculates the saturation vapor pressure over water for given temperature and pressure
#'@param temp temperature [K]
#'@return E_s, saturation vapor pressure over water [Pa]
#'@export
#'
calc_satvaporpressure = function(temp) {
    a=0.61094
    b=17.625
    c=243.04
    temp=temp-273.15 #K to deg C
    return(a*exp(b*temp/(temp+c))*1000)
}

#' Vapor pressure deficit (VPD)
#'
#'@description Calculates vapor pressure deficit (VPD) from temperature and relative humidity using Arrhenius formula
#'@param temp temperature [K]
#'@param rh relative humidity [percent]
#'@return VPD, vapor pressure deficit [Pa]
#'@export
#'
calc_vpd = function(temp,rh) {
    #calculate saturation pressure here using Arrhenius formula
    a=-1.0440397*10^4
    b=-11.29465
    c=-2.7022355*10^-2
    d=1.289036*10^-5
    e=-2.4780681*10^-9
    f=6.5459673
    temp=temp*5/9 #K to deg R (Rankine scale)
    es=exp(a/temp+b+c*temp+d*temp^2+e*temp^3+f*log(temp))
    return(es*(1-rh/100))
}

#' Potential temperature
#'
#'@description Calculates potential temperature for given temperature and pressure
#'@param temp temperature [K]
#'@param pres pressure [Pa]
#'@return potential temperature [K]
#'@export
#'
calc_theta = function(temp,pres) {
    return(temp*(100000/pres)^(Rd()/cp()))
}

#' Virtual temperature
#'
#'@description Calculates virtual temperature for given temperature and specific humidity (mixing ratio)
#'@param temp temperature [K]
#'@param q specific humidity [kg/kg]
#'@return virtual temperature [K]
#'@export
#'
calc_Tv = function(temp,q) {
    return(temp*(1+Rd()/Rv()*q))
}

#' Converts pressure to height (using barometric formula)
#'
#'@description Calculates height from pressure
#'@param pres pressure [Pa]
#'@param pres0 reference pressure, scalar [Pa], default \code{pres0=101315}
#'@param temp0 reference temperature, scalar [K], default \code{temp0=288.15}
#'@return height [m]
#'@export
#'
pres2height = function(pres,pres0=101315,temp0=288.15) {
    return(temp0/0.0065*(1-(pres/pres0)^(1/5.255)))
}

#' Converts relative humidity to specific humidity
#'
#'@description Calculates specific humidity from relative humidity, temperature and pressure
#'@param rh relative humidity [percent]
#'@param temp temperature [K]
#'@param pres pressure [Pa]
#'@return specific humidity [kg/kg]
#'@export
#'
rh2q = function(rh,temp,pres) {
    es=calc_satvaporpressure(temp) #saturation vapor pressure [Pa]
    e=rh*es/100 #vapor pressure [Pa]
    w=e*Rd()/(Rv()*(pres-e))
    return(w/(w+1))
}

#' Converts relative humidity to absolute humidity
#'
#'@description Calculates absolute humidity from relative humidity and temperature
#'@param rh relative humidity [percent]
#'@param temp temperature [K]
#'@return absolute humidity [kg/m^3]
#'@export
#'
rh2ah = function(rh,temp) {
    es=calc_satvaporpressure(temp)/100 #saturation vapor pressure [hPa, here]
    return(es*rh*2.1674/temp/1000)
}

#' Converts absolute humidity to relative humidity
#'
#'@description Calculates absolute humidity from relative humidity and temperature
#'@param ah absolute humidity [kg/m^3]
#'@param temp temperature [K]
#'@return relative humidity [percent]
#'@export
#'
ah2rh = function(ah,temp) {
    es=calc_satvaporpressure(temp)/100 #saturation vapor pressure [hPa, here]
    return(ah/(2.1674*es)*temp*1000)
}

#' Clear Sky Index (CSI)
#'
#'@description Calculates clear sky index
#'@param temp temperature [K]
#'@param lw_in longwave incoming radiation [W/m^2]
#'@param rh relative humidity [percent]
#'@param e vapor pressure [Pa] (either rh or e have to be given)
#'@return CSI, clear sky index
#'@export
#'
calc_csi = function(temp,lw_in,rh=NULL,e=NULL) {
    if (is.null(rh) & is.null(e)) {
        stop("Either relative humidity rh or vapor pressure e have to be given.")
    }
    if (!is.null(rh)) { #calculate vapor pressure
        es = calc_satvaporpressure(temp)
        e = rh * es/100
    }
    sigma=5.67*10^(-8)
    epsilon_A = lw_in/(sigma*temp^4) #actual atmospheric emissivity
    epsilon = 0.23 + 0.47*(100*e/temp)^(1/8) #(theoretical) clear sky emissivity, Marty and Philiponna, 2002
    return(epsilon_A/epsilon)
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
calc_windDirection = function(u,v) {
	return((180+180/pi*atan2(v,u))%%360) #from ERA5 doc: https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398
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
calc_gustfactor = function(ws_max,ws_mean) {
	return(ws_max/ws_mean)
}
