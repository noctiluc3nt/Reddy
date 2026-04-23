#' Eddy-covariance post-processing of virtual measurements from high-frequency model output
#'
#'@description Straightforward eddy-covariance post-processing routine for virtual measurements, assuming a given timeseries with the length of the desired averaging time
#'@param u u-wind [m/s]
#'@param v v-wind [m/s]
#'@param w w-wind [m/s]
#'@param theta potential temperature [K]
#'@param s scalar concentration [arbitrary unit] ( optional)
#'@param measurement_height measurement height [m], only used for calculation of the stability parameter \code{zeta}
#'@param do_detrending logical, should the data be linearly detrended? default \code{FALSE}
#'@param do_double_rotation logical, should the wind data be double rotated? default \code{TRUE}
#'@param do_flagging logical, should the data be flagged? default \code{TRUE}, i.e. several flags are calculated, but no data is removed, can be used for quality analysis
#'
#'@importFrom pracma detrend
#'
#'@return data frame of post-processed eddy-covariance data
#'@export
#'
#'
EC_processing_virtual_measurements = function(u,v,w,theta,s=NULL,measurement_height=1,do_detrending=FALSE,do_double_rotation=TRUE,do_flagging=TRUE) {
    #given data
    do_scalar=!(is.null(s))
    #check length of data, time resolution, desired averaging time and their consistence
    ndata=length(u)
    #prepare output data
    cat("\n... allocate storage for output data ...")
    column_names=c("u_mean","v_mean","w_mean","theta_mean","s_mean",
                    "u_sd","v_sd","w_sd","theta_sd","s_sd",
                    "wd","ws",
                    "tke","ustar","L","zeta","xb","yb",
                    "cov_uw","cov_vw","cov_uv","cov_thetaw","cov_sw",
                    "flag_all","flag_stationarity","flag_w","flag_itc",
                    "rotation_angle1","rotation_angle2",
                    "ampl_res_u","ampl_res_v","ampl_res_w","ampl_res_theta","ampl_res_s")
    out=array(NA,dim=c(1,length(column_names)))
    out=as.data.frame(out)
    colnames(out)=column_names
    #amplitude resolution
    out$ampl_res_u=get_amplitude_resolution(u)
    out$ampl_res_v=get_amplitude_resolution(v)
    out$ampl_res_w=get_amplitude_resolution(w)
    out$ampl_res_theta=get_amplitude_resolution(theta)
    if (do_scalar) out$ampl_res_s=get_amplitude_resolution(s)
    #wind (before rotation, assumes that the sonic is oriented towards north as indicated on the instrument)
    out$ws=mean(calc_windspeed(u,v),na.rm=TRUE)
    out$wd=calc_circular_mean(calc_windDirection(u,v))
    #double rotation
    if (do_double_rotation==TRUE) {
        wind_rotated=rotate_double(u,v,w)
        out$rotation_angle1=wind_rotated$theta
        out$rotation_angle2=wind_rotated$phi
        u=wind_rotated$u
        v=wind_rotated$v
        w=wind_rotated$w
    } 
    #flagging: stationarity
    if (do_flagging) out$flag_stationarity=flag_stationarity(theta,w,nsub=as.integer(ndata/6))
    #detrending
    if (do_detrending == TRUE) {
        u=pracma::detrend(u)
        v=pracma::detrend(v)
        w=pracma::detrend(w)
        theta=pracma::detrend(theta)
        if (do_scalar) s=pracma::detrend(s)
    }
    #averaging
    cat("\n... do time averaging ...")
    out$u_mean=mean(u)
    out$u_sd=sd(u)
    out$v_mean=mean(v)
    out$v_sd=sd(v)
    out$w_mean=mean(w)
    out$w_sd=sd(w)
    out$theta_mean=mean(theta)
    out$theta_sd=sd(theta)
    if (do_scalar) out$s_mean=mean(s)
    if (do_scalar) out$s_sd=sd(s)
    out$cov_uw=cov(u,w)
    out$cov_uv=cov(u,v)
    out$cov_vw=cov(v,w)
    out$cov_thetaw=cov(w,theta)
    if (do_scalar) out$cov_sw=cov(s,w)
    #calculate turbulence statistics
    out$tke=calc_tke(out$u_sd,out$v_sd,out$w_sd)
    out$ustar=calc_ustar(out$cov_uw,out$cov_vw)
    out$L=calc_L(out$ustar,out$theta_mean,out$cov_thetaw)
    out$zeta=calc_zeta(measurement_height,out$L)
    aniso=calc_anisotropy(out$u_sd,out$cov_uv,out$cov_uw,out$v_sd,out$cov_vw,out$w_sd)
    out$xb=aniso$xb
    out$yb=aniso$yb
    #flagging: the other flags
    if (do_flagging==TRUE) {
        out$flag_itc=flag_most(out$w_sd,out$ustar,out$zeta)
        out$flag_w=flag_w(out$w_mean)
        out$flag_all=pmax(out$flag_stationarity,out$flag_itc,out$flag_w,na.rm=TRUE)
    }
    cat("\n... done :) \n\n")
    return(out)
}
