#### Constants #####
#' sigma
#'
#' Stefan-Boltzmann constant [W/m^2/K^4]
#' @keywords internal
sigma=function(){
	return(5.67*10^(-8))
}

#' cp
#'
#' heat capacity, constant pressure [J/(kg*K)]
#' @keywords internal
cp=function(){
	return(1005.7)
}

#' g
#'
#' gravitational accelaration [m/s^2]
#' @keywords internal
g=function(){
	return(9.8062)
}

#' clight
#'
#' speed of light [m/s]
#' @keywords internal
clight=function(){
	return(299792458)
}

#' von Karman constant
#'
#' @keywords internal
karman=function() {
	return(0.4)
}

#' R
#'
#' ideal gas constant [J/(mol*K)]
#' @keywords internal
Runiversal = function() {
	return(8.31451)
}

#' rhoAir
#'
#' molar density of air [mol/m^3]
#' @keywords internal
rhoAir = function() {
	return(1.164)
}

#' Lv
#'
#' latent heat of vaporization [J/kg]
#' @keywords internal
Lv = function() {
	return(2264.705*1000)
}

#' M_H2O
#'
#' Molar mass of water
#' @keywords internal
M_H2O = function() {
	return(0.01802)
}

#' M_CO2
#'
#' Molar mass of carbon dioxide
#' @keywords internal
M_CO2 = function() {
	return(0.044)
}

#' M_CH4
#'
#' Molar mass of methane
#' @keywords internal
M_CH4 = function() {
	return(0.01604)
}