#' Fit Arrhenius-type temperature-respiration model
#'
#'@description Fit Arrhenius-type temperature-respiration model, e.g. applied to fit for nighttime (when Reco = NEE) and use this further to estimate GPP
#'@param Reco ecosystem respiration (during nighttime: Reco = NEE, CO2-flux), timeseries
#'@param temp temprature [K], timeseries corresponding to Reco
#'@param Rref reference respiration [J/(mol*K)], default \code{Rref=8.314}
#'@param Tref reference temperature [K], default \code{T_ref=283.15}
#'
#'@return summary of non-linear fit (with nls) for the activation energy using the Arrhenius type-temperature model
#'@export
#'
fit_respiration_temperature_relation=function(Reco,temp,Rref=8.134,Tref=283.15) {
    R0=median(Reco,na.rm=TRUE)
    fit = nls(Reco ~ R0*exp(-activation_energy/Rref*(1/temp-1/Tref)),
        data=list("Reco"=Reco,"temp"=temp),
        start=list("R0"=R0,"activation_energy"=60000))
    return(fit)
}