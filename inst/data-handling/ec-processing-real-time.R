#' Eddy-covariance post-processing for near-real-time analysis
#'
#'@description An example for an eddy-covariance post-processing routine utilizing the functions from ec_processing.R
#'@param u u-wind [m/s] (sonic)
#'@param v v-wind [m/s] (sonic)
#'@param w w-wind [m/s] (sonic)
#'@param temp temperature [K] (sonic)
#'@param h2o H2O mixing ratio (gas analyzer, optional)
#'@param co2 CO2 mixing ratio (gas analyzer, optional)
#'@param ch4 CH4 mixing ratio (gas analyzer, optional)
#'@param time_resolution time resolution of the measurements [s], default 20 Hz = 0.05 s
#'@param time_averaging desired time averaging for flux calculations [min], default 30 minutes
#'@param measurement_height measurement height [m], only used for calculation of the stability parameter \code{zeta}
#'@param do_despiking logical, should the data be despiked? default \code{TRUE}
#'@param despike_u vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_u=c(-15,15,10,2,8)}
#'@param despike_v vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_v=c(-15,15,10,2,8)}
#'@param despike_w vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_w=c(-4,4,10,2,8)}
#'@param despike_temp vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_temp=c(230,300,10,2,8)}
#'@param despike_h2o vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_h2o=c(0,12,10,2,8)}
#'@param despike_co2 vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_co2=c(0,12,10,2,8)}
#'@param despike_ch4 vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_ch4=c(0,12,10,2,8)}
#'@param do_detrending logical, should the data be linearly detrended? default \code{FALSE}
#'@param do_double_rotation logical, should the wind data be double rotated? default \code{TRUE}
#'@param do_flagging logical, should the data be flagged? default \code{TRUE}, i.e. several flags are calculated, but no data is removed, can be used for quality analysis
#'@param dir_blocked vector containing 2 elements: wind directions blocked through mast or tower, used in flow distortion flag only
#'@param do_SNDcorrection logical, should SND correction be applied to the buoyancy flux? default \code{TRUE}
#'@param do_WPLcorrection logical, should WPL correction be applied to gas fluxes? default \code{FALSE}
#'@param A constant used in SND correction, default \code{A=7/8} for CSAT3 sonic
#'@param B constant used in SND correction, default \code{A=7/8} for CSAT3 sonic
#'@param store logical, should the output be stored? default \code{TRUE}
#'@param format_out file format of the output, can be either \code{txt} or \code{rds} (for netcdf, see separate function), only used if \code{store=TRUE}
#'@param filename desired output filename, default \code{NULL}, the date and time of the run will be used to create a filename, only used if \code{store=TRUE}
#'@param meta logical, should meta data be stored? default \code{TRUE}
#'
#'@importFrom pracma detrend
#'
#'@return data frame of post-processed eddy-covariance data (that is also stored in the output file by default)
#'@export
#'
#'
EC_processing_realtime = function(u,v,w,temp,h2o=NULL,co2=NULL,ch4=NULL,
    time_resolution=0.05, #s
    time_averaging=30, #mins
    measurement_height=1, #m
    do_despiking=TRUE,despike_u=c(-15,15,10,2,8),despike_v=c(-15,15,10,2,8),despike_w=c(-4,4,10,2,8),despike_temp=c(230,300,10,2,8),despike_h2o=NULL,despike_co2=NULL,despike_ch4=NULL,
    do_detrending=FALSE,
    do_double_rotation=TRUE,
    do_flagging=TRUE, dir_blocked=c(0,0),
    do_SNDcorrection=TRUE,A=7/8,B=7/8,
    do_WPLcorrection=FALSE,
    store=TRUE,format_out="txt",filename=NULL,
    meta=FALSE
    ) {
    #given data
    do_h2o=!(is.null(h2o))
    do_co2=!(is.null(co2))
    do_ch4=!(is.null(ch4))
    #check length of data, time resolution, desired averaging time and their consistence
    ndata=length(u)
    meas_time=ndata*time_resolution #total measurement time [s]
    nint=meas_time/(time_averaging*60) #number of output values
    if (nint<1) stop("The measurement time is too short for the desired averaging time -- check time resolution and desired averaging time.")
    nint=round(nint) #needs to be integer value
    lint=max(time_averaging)*60/(time_resolution) #length of averaging interval, i.e. number of measurements to be averaged
    #prepare output data
    cat("\n... allocate storage for output data ...")
    column_names=c("u_mean","v_mean","w_mean","Ts_mean","h2o_mean","co2_mean","ch4_mean",
                    "u_sd","v_sd","w_sd","Ts_sd","h2o_sd","co2_sd","ch4_sd",
                    "wd","ws",
                    "tke","ustar","L","zeta",
                    "cov_uw","cov_vw","cov_uv","cov_wTs","cov_vTs","cov_h2ow","cov_co2w","cov_ch4w",
                    "cov_wT_snd",
                    "sh","lh","co2_flux","ch4_flux",
                    "flag_all","flag_stationarity","flag_distortion","flag_w","flag_itc",
                    "rotation_angle1","rotation_angle2",
                    "nr_spikes_u","nr_spikes_v","nr_spikes_w","nr_spikes_Ts","nr_spikes_h2o","nr_spikes_co2","nr_spikes_ch4",
                    "ampl_res_u","ampl_res_v","ampl_res_w","ampl_res_Ts","ampl_res_h2o","ampl_res_co2","ampl_res_ch4")
    out=array(NA,dim=c(nint,length(column_names)))
    out=as.data.frame(out)
    colnames(out)=column_names
    #despiking
    if (do_despiking==TRUE) {
        cat("\n... do despiking ...")
        #todo add spike counting
        u=despiking(u,c(despike_u[1],despike_u[2]),despike_u[3],despike_u[4],despike_u[5])
        v=despiking(v,c(despike_v[1],despike_v[2]),despike_v[3],despike_v[4],despike_v[5])
        w=despiking(w,c(despike_w[1],despike_w[2]),despike_w[3],despike_w[4],despike_w[5])
        temp=despiking(temp,c(despike_temp[1],despike_temp[2]),despike_temp[3],despike_temp[4],despike_temp[5])
        out$nr_spikes_u=count_spikes(u,c(despike_u[1],despike_u[2]))
        out$nr_spikes_v=count_spikes(v,c(despike_v[1],despike_v[2]))
        out$nr_spikes_w=count_spikes(w,c(despike_w[1],despike_w[2]))
        out$nr_spikes_Ts=count_spikes(temp,c(despike_temp[1],despike_temp[2]))
        if (do_h2o & !is.null(despike_h2o)) h2o=despiking(h2o,c(despike_h2o[1],despike_h2o[2]),despike_h2o[3],despike_h2o[4],despike_h2o[5])
        if (do_co2 & !is.null(despike_co2)) co2=despiking(co2,c(despike_co2[1],despike_co2[2]),despike_co2[3],despike_co2[4],despike_co2[5])
        if (do_ch4 & !is.null(despike_ch4)) ch4=despiking(ch4,c(despike_ch4[1],despike_ch4[2]),despike_ch4[3],despike_ch4[4],despike_ch4[5])
        if (do_h2o & !is.null(despike_h2o)) out$nr_spikes_h2o=count_spikes(h2o,c(despike_h2o[1],despike_h2o[2]))
        if (do_co2 & !is.null(despike_co2)) out$nr_spikes_co2=count_spikes(co2,c(despike_co2[1],despike_co2[2]))
        if (do_ch4 & !is.null(despike_ch4)) out$nr_spikes_ch4=count_spikes(ch4,c(despike_ch4[1],despike_ch4[2]))
    }
    #amplitude resolution
    out$ampl_res_u=get_amplitude_resolution(u)
    out$ampl_res_v=get_amplitude_resolution(v)
    out$ampl_res_w=get_amplitude_resolution(w)
    out$ampl_res_Ts=get_amplitude_resolution(temp)
    if (do_h2o) out$ampl_res_h2o=get_amplitude_resolution(h2o)
    if (do_co2) out$ampl_res_co2=get_amplitude_resolution(co2)
    if (do_ch4) out$ampl_res_ch4=get_amplitude_resolution(ch4)
    cat("\n... start loop over data: do double rotation and stationarity flagging (if requested) ...")
    for (i in 1:nint) {
        #cat(paste0("\n\t #index: ",i,"\t progress: ",round(i/nint*100,2)," %"))
        i1=(i-1)*lint+1
        i2=(i*lint)
        iselect=seq(i1,i2)
        #wind (before rotation, assumes that the sonic is oriented towards north as indicated on the instrument)
        out$ws[i]=mean(calc_windspeed(u[iselect],v[iselect]),na.rm=TRUE)
        out$wd[i]=calc_circular_mean(calc_windDirection(u[iselect],v[iselect]))
        if (do_flagging) out$flag_distortion[i]=flag_distortion(u[iselect],v[iselect],dir_blocked)
        #double rotation
        if (do_double_rotation==TRUE) {
            wind_rotated=rotate_double(u[iselect],v[iselect],w[iselect])
            out$rotation_angle1[i]=wind_rotated$theta
            out$rotation_angle2[i]=wind_rotated$phi
            u[iselect]=wind_rotated$u
            v[iselect]=wind_rotated$v
            w[iselect]=wind_rotated$w
        } 
        #flagging: stationarity
        if (do_flagging) out$flag_stationarity[i]=flag_stationarity(temp[iselect],w[iselect],nsub=as.integer(lint/4))
    }
    #detrending
    if (do_detrending == TRUE) {
        u=pracma::detrend(u)
        v=pracma::detrend(v)
        w=pracma::detrend(w)
        temp=pracma::detrend(temp)
        if (do_h2o) h2o=pracma::detrend(h2o)
        if (do_co2) co2=pracma::detrend(co2)
        if (do_ch4) ch4=pracma::detrend(ch4)
    }
    #averaging
    cat("\n... do time averaging ...")
    u_avg=accumulate_timeseries(u,time_resolution,time_averaging*60)
    out$u_mean=u_avg$mean[[1]]
    out$u_sd=u_avg$sd[[1]]
    v_avg=accumulate_timeseries(v,time_resolution,time_averaging*60)
    out$v_mean=v_avg$mean[[1]]
    out$v_sd=v_avg$sd[[1]]
    w_avg=accumulate_timeseries(w,time_resolution,time_averaging*60)
    out$w_mean=w_avg$mean[[1]]
    out$w_sd=w_avg$sd[[1]]
    Ts_avg=accumulate_timeseries(temp,time_resolution,time_averaging*60)
    out$Ts_mean=Ts_avg$mean[[1]]
    out$Ts_sd=Ts_avg$sd[[1]]
    if (do_h2o) {
        h2o_avg=accumulate_timeseries(h2o,time_resolution,time_averaging*60)
        out$h2o_mean=h2o_avg$mean[[1]]
        out$h2o_sd=h2o_avg$sd[[1]]
    }
    if (do_co2) {
        co2_avg=accumulate_timeseries(co2,time_resolution,time_averaging*60)
        out$co2_mean=co2_avg$mean[[1]]
        out$co2_sd=co2_avg$sd[[1]]
    }
    if (do_ch4) {
        ch4_avg=accumulate_timeseries(ch4,time_resolution,time_averaging*60)
        out$ch4_mean=ch4_avg$mean[[1]]
        out$ch4_sd=ch4_avg$sd[[1]]
    }
    #flux calculation
    cat("\n... do flux calculation ...")
    out$cov_uw=accumulate_timeseries(u*w,time_resolution,time_averaging*60)$mean[[1]]-out$u_mean*out$w_mean
    out$cov_uv=accumulate_timeseries(u*v,time_resolution,time_averaging*60)$mean[[1]]-out$u_mean*out$v_mean
    out$cov_vw=accumulate_timeseries(v*w,time_resolution,time_averaging*60)$mean[[1]]-out$v_mean*out$w_mean
    out$cov_wTs=accumulate_timeseries(w*temp,time_resolution,time_averaging*60)$mean[[1]]-out$Ts_mean*out$w_mean
    out$cov_vTs=accumulate_timeseries(v*temp,time_resolution,time_averaging*60)$mean[[1]]-out$Ts_mean*out$v_mean
    if (do_h2o) out$cov_h2ow=accumulate_timeseries(w*h2o,time_resolution,time_averaging*60)$mean[[1]]-out$h2o_mean*out$w_mean
    if (do_co2) out$cov_co2w=accumulate_timeseries(w*co2,time_resolution,time_averaging*60)$mean[[1]]-out$co2_mean*out$w_mean
    if (do_ch4) out$cov_ch4w=accumulate_timeseries(w*ch4,time_resolution,time_averaging*60)$mean[[1]]-out$ch4_mean*out$w_mean
    #SND correction
    if (do_SNDcorrection==TRUE) {
        if (do_h2o == FALSE) { #cross-wind correction only
            out$cov_wT_snd=SNDcorrection(out$Ts_mean,out$u_mean,out$v_mean,out$cov_uw,out$cov_vw,out$cov_wTs,NULL,A,B)
        } else {
            out$cov_wT_snd=SNDcorrection(out$Ts_mean,out$u_mean,out$v_mean,out$cov_uw,out$cov_vw,out$cov_wTs,out$cov_h2ow,A,B)
        }
        out$sh=cov2sh(out$cov_wT_snd)
    } else {
        out$sh=cov2sh(out$cov_wTs)
    }
    if (do_h2o) out$lh=cov2lh(out$cov_h2ow)
    if (do_co2) out$co2_flux=cov2cf(out$cov_co2w)
    #WPL correction
    if (do_WPLcorrection==TRUE) {
        if (do_h2o) out$cov_h2ow=WPLcorrection(out$Ts_mean,out$h2o_mean/1.2,out$cov_wTs,out$h2o_mean,out$cov_h2ow)
        if (do_co2) out$cov_co2w=WPLcorrection(out$Ts_mean,out$h2o_mean/1.2,out$cov_wTs,out$h2o_mean,out$cov_h2ow,out$co2_mean,out$cov_co2w)
        if (do_ch4) out$cov_ch4w=WPLcorrection(out$Ts_mean,out$h2o_mean/1.2,out$cov_wTs,out$h2o_mean,out$cov_h2ow,out$ch4_mean,out$cov_ch4w)
        out$lh=cov2lh(out$cov_h2ow)
        out$co2_flux=cov2cf(out$cov_co2w)
    }
    #calculate turbulence statistics
    out$tke=calc_tke(out$u_sd,out$v_sd,out$w_sd)
    out$ustar=calc_ustar(out$cov_uw,out$cov_vw)
    out$L=calc_L(out$ustar,out$Ts_mean,out$cov_wTs)
    out$zeta=calc_zeta(measurement_height,out$L)
    #flagging: the other flags
    if (do_flagging==TRUE) {
        out$flag_itc=flag_most(out$w_sd,out$ustar,out$zeta)
        out$flag_w=flag_w(out$w_mean)
        out$flag_all=pmax(out$flag_stationarity,out$flag_distortion,out$flag_itc,out$flag_w,na.rm=TRUE)
    }
    #------------------------------------------------
    #store post-processed data
    systime=Sys.time()
    systime_string=format(systime,"%F_%H%M%S",tz="utc")
    if (is.null(filename)) {
        filename=paste0("ec-processing_Reddy_",systime_string,".",format_out)
    }
    #out=out[,colSums(is.na(out))<nint] #remove columns that contain only NA
    if (store == TRUE) {    
        if (format_out=="txt" | format_out=="dat") {
            cat("\n... store output as .dat file ...")
            write.table(out,file=filename,quote=FALSE,row.names=FALSE,col.names=TRUE)
        } else if (format_out=="csv") {
            cat("\n... store output as .csv file ...")
            write.csv(out,file=filename,quote=FALSE,row.names=FALSE,col.names=TRUE)
        } else if (format_out=="rds") {
            cat("\n... store output as .rds file ...")
            saveRDS(out,file=filename)
        #} else if (format_out=="nc" | format_out=="netcdf") {
        #   #see separate function in data_handling directory: the storage as netcdf is not available here, since installing the ncdf4 package sometimes causes problems, but the respective function ECprocessing_nc can be found in the directory data_handling
        } else {
            warning("You chose an output format that is not available by default in this function. But you can store the returned data frame in the desired format yourself.")
        }
    }
    #write metadata file
    if (meta) {
        meta=paste0("--------------------------------------------------\nEddy-covariance post-processing with Reddy package\n--------------------------------------------------
             \ndate: ",format(systime,"%F %T",tz="utc")," UTC",
            "\noutput filename: ", filename,
            "\ntime resolution of input data: ", time_resolution, " s",
            "\naveraging time: ", time_averaging, " min",
            "\ndo_h2o: ", do_h2o,
            "\ndo_co2: ", do_co2,
            "\ndo_ch4: ", do_ch4,
            "\ndo_despiking: ", do_despiking,
            "\ndo_detrending: ", do_detrending,
            "\ndo_double_rotation: ", do_double_rotation,
            #"\ndo_planar_fit: ", do_planar_fit,
            "\ndo_flagging: ", do_flagging,
            "\ndo_SNDcorrection: ", do_SNDcorrection
        )
        cat(paste0("\n\n",meta,"\n\n"))
        if (store==TRUE) writeLines(meta,paste0("metadata_ec-processing_Reddy_",systime_string,".txt"))
    }
    return(out)
}


#' Apply quality control on high-frequency data (e.g. as post-processing for MRD or quadrant analysis)
#'
#'@description Applies quality control and rotation to high-frequency data (and outputs the high-frequency data)
#'@param u u-wind [m/s] (sonic)
#'@param v v-wind [m/s] (sonic)
#'@param w w-wind [m/s] (sonic)
#'@param temp temperature [K] (sonic)
#'@param h2o H2O mixing ratio (gas analyzer, optional)
#'@param co2 CO2 mixing ratio (gas analyzer, optional)
#'@param ch4 CH4 mixing ratio (gas analyzer, optional)
#'@param do_despiking logical, should the data be despiked? default \code{TRUE}
#'@param despike_u vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_u=c(-15,15,10,2,8)}
#'@param despike_v vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_v=c(-15,15,10,2,8)}
#'@param despike_w vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_w=c(-4,4,10,2,8)}
#'@param despike_temp vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_temp=c(230,300,10,2,8)}
#'@param despike_h2o vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_h2o=c(0,12,10,2,8)}
#'@param despike_co2 vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_co2=c(0,12,10,2,8)}
#'@param despike_ch4 vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_ch4=c(0,12,10,2,8)}
#'@param do_double_rotation logical, should the wind data be double rotated? default \code{TRUE}
#'@param do_planar_fit logical, should the data be rotated with planar fit? default \code{FALSE} (either double rotation or planar fit can be \code{TRUE})
#'@param do_detrending logical, should the data be linearly detrended? default \code{FALSE}
#'@param do_flagging logical, should the data be flagged? default \code{TRUE}, i.e. several flags are calculated, but no data is removed, can be used for quality analysis
#'
#'@importFrom pracma detrend
#'
#'@return quality checked data in the same dimensions as the input variables
#'@export
#'
#'
apply_quality_control = function(u,v,w,temp,h2o=NULL,co2=NULL,ch4=NULL,
    do_despiking=TRUE,despike_u=c(-15,15,10,2,8),despike_v=c(-15,15,10,2,8),despike_w=c(-4,4,10,2,8),despike_temp=c(230,300,10,2,8),despike_h2o=c(0,12,10,2,8),despike_co2=c(300,500,10,4,10),despike_ch4=c(0,12,10,2,8),
    do_double_rotation=TRUE,
    do_planar_fit=FALSE,
    do_detrending=FALSE,
    do_flagging=TRUE
    ) {
    #given data
    do_h2o=!(is.null(h2o))
    do_co2=!(is.null(co2))
    do_ch4=!(is.null(ch4))
    ndata=length(u)
    #despiking
    if (do_despiking==TRUE) {
        cat("\n... do despiking ...")
        u=despiking(u,c(despike_u[1],despike_u[2]),despike_u[3],despike_u[4],despike_u[5])
        v=despiking(v,c(despike_v[1],despike_v[2]),despike_v[3],despike_v[4],despike_v[5])
        w=despiking(w,c(despike_w[1],despike_w[2]),despike_w[3],despike_w[4],despike_w[5])
        temp=despiking(temp,c(despike_temp[1],despike_temp[2]),despike_temp[3],despike_temp[4],despike_temp[5])
        if (do_h2o) h2o=despiking(h2o,c(despike_h2o[1],despike_h2o[2]),despike_h2o[3],despike_h2o[4],despike_h2o[5])
        if (do_co2) co2=despiking(co2,c(despike_co2[1],despike_co2[2]),despike_co2[3],despike_co2[4],despike_co2[5])
        if (do_ch4) ch4=despiking(ch4,c(despike_ch4[1],despike_ch4[2]),despike_ch4[3],despike_ch4[4],despike_ch4[5])
    }
    #rotation
    if (do_double_rotation == TRUE) {
        wind_rotated=rotate_double(u,v,w)
        u=wind_rotated$u
        v=wind_rotated$v
        w=wind_rotated$w
    }
    if (do_planar_fit == TRUE) {
        wind_rotated=rotate_planar(u,v,w)
        u=wind_rotated$u
        v=wind_rotated$v
        w=wind_rotated$w    
    }
    #detrending
    if (do_detrending == TRUE) {
        u=pracma::detrend(u)
        v=pracma::detrend(v)
        w=pracma::detrend(w)
        temp=pracma::detrend(temp)
        if (do_h2o) h2o=pracma::detrend(h2o)
        if (do_co2) co2=pracma::detrend(co2)
        if (do_ch4) ch4=pracma::detrend(ch4)
    }
    #------------------------------------------------
    # prepare return 
    out=list(u=u,v=v,w=w,Ts=temp)
    if(do_h2o) out$h2o=h2o
    if(do_co2) out$co2=co2
    if(do_ch4) out$ch4=ch4
    return(out)
}


########################################################3

#' Eddy-covariance post-processing for near-real-time analysis
#'
#'@description An example for an eddy-covariance post-processing routine utilizing the functions from ec_processing.R
#'@param u u-wind [m/s] (sonic)
#'@param v v-wind [m/s] (sonic)
#'@param w w-wind [m/s] (sonic)
#'@param temp temperature [K] (sonic)
#'@param h2o H2O mixing ratio (gas analyzer, optional)
#'@param co2 CO2 mixing ratio (gas analyzer, optional)
#'@param ch4 CH4 mixing ratio (gas analyzer, optional)
#'@param measurement_height measurement height [m], only used for calculation of the stability parameter \code{zeta}
#'@param do_despiking logical, should the data be despiked? default \code{TRUE}
#'@param despike_u vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_u=c(-15,15,10,2,8)}
#'@param despike_v vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_v=c(-15,15,10,2,8)}
#'@param despike_w vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_w=c(-4,4,10,2,8)}
#'@param despike_temp vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_temp=c(230,300,10,2,8)}
#'@param despike_h2o vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_h2o=c(0,12,10,2,8)}
#'@param despike_co2 vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_co2=c(0,12,10,2,8)}
#'@param despike_ch4 vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_ch4=c(0,12,10,2,8)}
#'@param do_detrending logical, should the data be linearly detrended? default \code{FALSE}
#'@param do_double_rotation logical, should the wind data be double rotated? default \code{TRUE}
#'@param do_flagging logical, should the data be flagged? default \code{TRUE}, i.e. several flags are calculated, but no data is removed, can be used for quality analysis
#'@param dir_blocked vector containing 2 elements: wind directions blocked through mast or tower, used in flow distortion flag only
#'@param do_SNDcorrection logical, should SND correction be applied to the buoyancy flux? default \code{TRUE}
#'@param do_WPLcorrection logical, should WPL correction be applied to gas fluxes? default \code{FALSE}
#'@param do_shift2maxccf logical, should the time series (per averaging interval) be shifted to maximize their cross-correlation? default \code{TRUE}
#'@param maxlag maximum lag used in \code{shift2maxccf()}, default 20 (i.e. 1 s when sampling in 20 Hz)
#'@param A constant used in SND correction, default \code{A=7/8} for CSAT3 sonic
#'@param B constant used in SND correction, default \code{A=7/8} for CSAT3 sonic
#'
#'@importFrom pracma detrend
#'
#'@return data frame of post-processed eddy-covariance data (that is also stored in the output file by default)
#'@export
#'
#'
EC_processing_per_file = function(u,v,w,temp,h2o=NULL,co2=NULL,ch4=NULL,
    measurement_height=1, #m
    do_despiking=TRUE,despike_u=c(-15,15,10,2,8),despike_v=c(-15,15,10,2,8),despike_w=c(-4,4,10,2,8),despike_temp=c(230,300,10,2,8),despike_h2o=NULL,despike_co2=NULL,despike_ch4=NULL,
    do_detrending=FALSE,
    do_double_rotation=TRUE,
    do_flagging=TRUE, dir_blocked=c(0,0),
    do_SNDcorrection=TRUE,A=7/8,B=7/8,
    do_WPLcorrection=FALSE,
    do_shift2maxccf=TRUE,maxlag=20,
    store=TRUE,format_out="txt",filename=NULL,
    meta=FALSE
    ) {
    #given data
    do_h2o=!(is.null(h2o))
    do_co2=!(is.null(co2))
    do_ch4=!(is.null(ch4))
    #prepare output data
    cat("\n... allocate storage for output data ...")
    column_names=c("u_mean","v_mean","w_mean","Ts_mean","h2o_mean","co2_mean","ch4_mean",
                    "u_sd","v_sd","w_sd","Ts_sd","h2o_sd","co2_sd","ch4_sd",
                    "wd","ws",
                    "tke","ustar","L","zeta",
                    "cov_uw","cov_vw","cov_uv","cov_wTs","cov_vTs","cov_h2ow","cov_co2w","cov_ch4w",
                    "cov_wT_snd",
                    "sh","lh","co2_flux","ch4_flux",
                    "flag_all","flag_stationarity","flag_distortion","flag_w","flag_itc",
                    "rotation_angle1","rotation_angle2",
                    "nr_spikes_u","nr_spikes_v","nr_spikes_w","nr_spikes_Ts","nr_spikes_h2o","nr_spikes_co2","nr_spikes_ch4",
                    "ampl_res_u","ampl_res_v","ampl_res_w","ampl_res_Ts","ampl_res_h2o","ampl_res_co2","ampl_res_ch4")
    out=array(NA,dim=c(1,length(column_names)))
    out=as.data.frame(out)
    colnames(out)=column_names
    #despiking
    if (do_despiking==TRUE) {
        cat("\n... do despiking ...")
        #todo add spike counting
        u=despiking(u,c(despike_u[1],despike_u[2]),despike_u[3],despike_u[4],despike_u[5])
        v=despiking(v,c(despike_v[1],despike_v[2]),despike_v[3],despike_v[4],despike_v[5])
        w=despiking(w,c(despike_w[1],despike_w[2]),despike_w[3],despike_w[4],despike_w[5])
        temp=despiking(temp,c(despike_temp[1],despike_temp[2]),despike_temp[3],despike_temp[4],despike_temp[5])
        out$nr_spikes_u=count_spikes(u,c(despike_u[1],despike_u[2]))
        out$nr_spikes_v=count_spikes(v,c(despike_v[1],despike_v[2]))
        out$nr_spikes_w=count_spikes(w,c(despike_w[1],despike_w[2]))
        out$nr_spikes_Ts=count_spikes(temp,c(despike_temp[1],despike_temp[2]))
        if (do_h2o & !is.null(despike_h2o)) h2o=despiking(h2o,c(despike_h2o[1],despike_h2o[2]),despike_h2o[3],despike_h2o[4],despike_h2o[5])
        if (do_co2 & !is.null(despike_co2)) co2=despiking(co2,c(despike_co2[1],despike_co2[2]),despike_co2[3],despike_co2[4],despike_co2[5])
        if (do_ch4 & !is.null(despike_ch4)) ch4=despiking(ch4,c(despike_ch4[1],despike_ch4[2]),despike_ch4[3],despike_ch4[4],despike_ch4[5])
        if (do_h2o & !is.null(despike_h2o)) out$nr_spikes_h2o=count_spikes(h2o,c(despike_h2o[1],despike_h2o[2]))
        if (do_co2 & !is.null(despike_co2)) out$nr_spikes_co2=count_spikes(co2,c(despike_co2[1],despike_co2[2]))
        if (do_ch4 & !is.null(despike_ch4)) out$nr_spikes_ch4=count_spikes(ch4,c(despike_ch4[1],despike_ch4[2]))
    }
    #amplitude resolution
    out$ampl_res_u=get_amplitude_resolution(u)
    out$ampl_res_v=get_amplitude_resolution(v)
    out$ampl_res_w=get_amplitude_resolution(w)
    out$ampl_res_Ts=get_amplitude_resolution(temp)
    if (do_h2o) out$ampl_res_h2o=get_amplitude_resolution(h2o)
    if (do_co2) out$ampl_res_co2=get_amplitude_resolution(co2)
    if (do_ch4) out$ampl_res_ch4=get_amplitude_resolution(ch4)
    #wind (before rotation, assumes that the sonic is oriented towards north as indicated on the instrument)
    out$ws=mean(calc_windspeed(u,v),na.rm=TRUE)
    out$wd=calc_circular_mean(calc_windDirection(u,v))
    if (do_flagging) out$flag_distortion=flag_distortion(u,v,dir_blocked)
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
    #if (do_flagging) out$flag_stationarity=flag_stationarity(temp,w)
    #detrending
    if (do_detrending == TRUE) {
        u=pracma::detrend(u)
        v=pracma::detrend(v)
        w=pracma::detrend(w)
        temp=pracma::detrend(temp)
        if (do_h2o) h2o=pracma::detrend(h2o)
        if (do_co2) co2=pracma::detrend(co2)
        if (do_ch4) ch4=pracma::detrend(ch4)
    }
    #averaging
    cat("\n... do time averaging ...")
    out$u_mean=mean(u,na.rm=T)
    out$u_sd=sd(u,na.rm=T)
    out$v_mean=mean(v,na.rm=T)
    out$v_sd=sd(v,na.rm=T)
    out$w_mean=mean(w,na.rm=T)
    out$w_sd=sd(w,na.rm=T)
    out$Ts_mean=mean(temp,na.rm=T)
    out$Ts_sd=sd(temp,na.rm=T)
    if (do_h2o) {
        out$h2o_mean=mean(h2o,na.rm=T)
        out$h2o_sd=sd(h2o,na.rm=T)
    }
    if (do_co2) {
        out$co2_mean=mean(co2,na.rm=T)
        out$co2_sd=sd(co2,na.rm=T)
    }
    if (do_ch4) {
        out$ch4_mean=mean(ch4,na.rm=T)
        out$ch4_sd=sd(ch4,na.rm=T)
    }
    #flux calculation
    cat("\n... do flux calculation ...")
    if (do_shift2maxccf==TRUE) {
        shifted=shift2maxccf(u,w,maxlag=maxlag)
        out$cov_uw=cov(shifted[,1],shifted[,2],use="pairwise.complete.obs")
        shifted=shift2maxccf(u,v,maxlag=maxlag)
        out$cov_uv=cov(shifted[,1],shifted[,2],use="pairwise.complete.obs")
        shifted=shift2maxccf(v,w,maxlag=maxlag)
        out$cov_vw=cov(shifted[,1],shifted[,2],use="pairwise.complete.obs")
        shifted=shift2maxccf(temp,w,maxlag=maxlag)
        out$cov_wTs=cov(shifted[,1],shifted[,2],use="pairwise.complete.obs")
        shifted=shift2maxccf(temp,v,maxlag=maxlag)
        out$cov_vTs=cov(shifted[,1],shifted[,2],use="pairwise.complete.obs")
        if (do_h2o) {
            shifted=shift2maxccf(h2o,w,maxlag=maxlag)
            out$cov_co2w=cov(shifted[,1],shifted[,2],use="pairwise.complete.obs")
        }
        if (do_co2) {
            shifted=shift2maxccf(co2,w,maxlag=maxlag)
            out$cov_co2w=cov(shifted[,1],shifted[,2],use="pairwise.complete.obs")
        }
        if (do_ch4) {
            shifted=shift2maxccf(ch4,w,maxlag=maxlag)
            out$cov_ch4w=cov(shifted[,1],shifted[,2],use="pairwise.complete.obs")
        }
    } else {
        out$cov_uw=cov(u,w,use="pairwise.complete.obs")
        out$cov_uv=cov(u,v,use="pairwise.complete.obs")
        out$cov_vw=cov(v,w,use="pairwise.complete.obs")
        out$cov_wTs=cov(w,temp,use="pairwise.complete.obs")
        out$cov_vTs=cov(v,temp,use="pairwise.complete.obs")
        if (do_h2o) {
            out$cov_h2ow=cov(h2o,w,use="pairwise.complete.obs")
        }
        if (do_co2) {
            out$cov_co2w=cov(co2,w,use="pairwise.complete.obs")
        }
        if (do_ch4) {
            out$cov_ch4w=cov(ch4,w,use="pairwise.complete.obs")
        }
    }
    #SND correction
    if (do_SNDcorrection==TRUE) {
        if (do_h2o == FALSE) { #cross-wind correction only
            out$cov_wT_snd=SNDcorrection(out$Ts_mean,out$u_mean,out$v_mean,out$cov_uw,out$cov_vw,out$cov_wTs,NULL,A,B)
        } else {
            out$cov_wT_snd=SNDcorrection(out$Ts_mean,out$u_mean,out$v_mean,out$cov_uw,out$cov_vw,out$cov_wTs,out$cov_h2ow,A,B)
        }
        out$sh=cov2sh(out$cov_wT_snd)
    } else {
        out$sh=cov2sh(out$cov_wTs)
    }
    if (do_h2o) out$lh=cov2lh(out$cov_h2ow)
    if (do_co2) out$co2_flux=cov2cf(out$cov_co2w)
    #WPL correction
    if (do_WPLcorrection==TRUE) {
        if (do_h2o) out$cov_h2ow=WPLcorrection(out$Ts_mean,out$h2o_mean/1.2,out$cov_wTs,out$h2o_mean,out$cov_h2ow)
        if (do_co2) out$cov_co2w=WPLcorrection(out$Ts_mean,out$h2o_mean/1.2,out$cov_wTs,out$h2o_mean,out$cov_h2ow,out$co2_mean,out$cov_co2w)
        if (do_ch4) out$cov_ch4w=WPLcorrection(out$Ts_mean,out$h2o_mean/1.2,out$cov_wTs,out$h2o_mean,out$cov_h2ow,out$ch4_mean,out$cov_ch4w)
        out$lh=cov2lh(out$cov_h2ow)
        out$co2_flux=cov2cf(out$cov_co2w)
    }
    #calculate turbulence statistics
    out$tke=calc_tke(out$u_sd,out$v_sd,out$w_sd)
    out$ustar=calc_ustar(out$cov_uw,out$cov_vw)
    out$L=calc_L(out$ustar,out$Ts_mean,out$cov_wTs)
    out$zeta=calc_zeta(measurement_height,out$L)
    #flagging: the other flags
    if (do_flagging==TRUE) {
        out$flag_itc=flag_most(out$w_sd,out$ustar,out$zeta)
        out$flag_w=flag_w(out$w_mean)
        out$flag_all=pmax(out$flag_stationarity,out$flag_distortion,out$flag_itc,out$flag_w,na.rm=TRUE)
    }
    #------------------------------------------------
    #store post-processed data
    systime=Sys.time()
    systime_string=format(systime,"%F_%H%M%S",tz="utc")
    if (is.null(filename)) {
        filename=paste0("ec-processing_Reddy_",systime_string,".",format_out)
    }
    return(out)
}

