#' Fit Arrhenius-type temperature-respiration model
#'
#'@description Fit Arrhenius-type temperature-respiration model, e.g. applied to fit for nighttime (when Reco = NEE) and use this further to estimate GPP
#'@param Reco ecosystem respiration [mu mol CO2/(m^2 s)] (during nighttime: Reco = NEE, CO2-flux), timeseries
#'@param temp temperature [K], timeseries corresponding to Reco
#'@param Rref reference respiration [J/(mol*K)], default \code{Rref=8.314}
#'@param Tref reference temperature [K], default \code{T_ref=283.15}
#'
#'@return summary of non-linear fit (with nls) for the activation energy using the Arrhenius type-temperature model
#'@export
#'
fit_respiration_temperature_relation=function(Reco,temp,Rref=8.134,Tref=283.15) {
    if (length(Reco) != length(temp)) {
        stop("The timeseries of Reco and temp have to be of the same length.")
    }
    Reco0=median(Reco,na.rm=TRUE)
    fit = nls(Reco ~ R0*exp(-activation_energy/Rref*(1/temp-1/Tref)),
        data=list("Reco"=Reco,"temp"=temp),
        start=list("R0"=Reco0,"activation_energy"=60000))
    return(fit)
}


#' Fit light-response curve
#'
#'@description Fit light-response curve (rectangular-hyperbola) of NEE based on PAR (for daytime conditions)
#'@param NEE net ecosystem respiration [mu mol CO2/(m^2 s)] (during daytime), timeseries
#'@param PAR PAR (or shortwave incoming radiation) [mu mol photons/(m^2 s)], timeseries corresponding to NEE
#'
#'@return summary of non-linear fit (with nls) for the activation energy using the Arrhenius type-temperature model
#'@export
#'
fit_light_response_curve=function(NEE,PAR) {
    if (length(NEE) != length(PAR)) {
        stop("The timeseries of NEE and PAR have to be of the same length.")
    }
    Reco0=median(NEE,na.rm=TRUE)
    alpha0=0.01
    Amax0=abs(min(NEE,na.rm=TRUE))
    fit=nls(NEE ~ Reco-(alpha*PAR*Amax)/(alpha*PAR+Amax),
        data=list("NEE"=NEE,"PAR"=PAR),
        start = list("Reco"=Reco0,"alpha"=alpha0,"Amax"=Amax0))
    return(fit)
}