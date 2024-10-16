#' Eddy-covariance post-processing
#'
#'@description An eddy-covariance post-processing routine utilizing the functions from ec_processing.R
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
#'@param do_despiking locigal, should the data be despiked? default \code{TRUE}
#'@param despike_u vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_u=c(-15,15,10,2,8)}
#'@param despike_v vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_v=c(-15,15,10,2,8)}
#'@param despike_w vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_w=c(-4,4,10,2,8)}
#'@param despike_temp vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_temp=c(230,300,10,2,8)}
#'@param despike_h2o vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_h2o=c(0,12,10,2,8)}
#'@param despike_co2 vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_co2=c(0,12,10,2,8)}
#'@param despike_ch4 vector containing 5 elements: lower and upper bound, MAD factor, threshold skewness, threshold kurtosis. Details see \code{?despiking}. Default \code{despike_ch4=c(0,12,10,2,8)}
#'@param do_detrending logical, should the data be linearly detrended? default \code{FALSE}
#'@param do_double_rotation locigal, should the wind data be double rotated? default \code{TRUE}
#'@param do_planar_fit locigal, should the data be rotated with planar fit? default \code{FALSE} (either double rotation or planar fit can be \code{TRUE})
#'@param do_flagging locigal, should the data be flagged? default \code{TRUE}, i.e. several flags are calculated, but no data is removed, can be used for quality analysis
#'@param dir_blocked vector containing 2 elements: wind directions blocked through mast or tower, used in flow distortion flag only
#'@param do_SNDcorrection locigal, should SND correction be applied to the buoyancy flux? default \code{TRUE}
#'@param A constant used in SND correction, default \code{A=7/8} for CSAT3 sonic
#'@param B constant used in SND correction, default \code{A=7/8} for CSAT3 sonic
#'@param store logical, should the output be stored? default \code{TRUE}
#'@param format_out file format of the output, can be either \code{txt} or \code{rds} (for netcdf, see separate function), only used if \code{store=TRUE}
#'@param filename desired output filename, default \code{NULL}, the date and runtime will be used to create a filename, only used if \code{store=TRUE}
#'@param meta locical, should meta data be stored? default \code{TRUE}
#'
#'@importFrom pracma detrend
#'
#'@return data frame of post-processed eddy-covariance data (that is also stored in the output file by default)
#'@export
#'
#'
ECprocessing = function(u,v,w,temp,h2o=NULL,co2=NULL,ch4=NULL,
    time_resolution=0.05, #s
    time_averaging=30, #mins
    measurement_height=1, #m
    do_despiking=TRUE,despike_u=c(-15,15,10,2,8),despike_v=c(-15,15,10,2,8),despike_w=c(-4,4,10,2,8),despike_temp=c(230,300,10,2,8),despike_h2o=c(0,12,10,2,8),despike_co2=c(0,12,10,2,8),despike_ch4=c(0,12,10,2,8),
    do_detrending=FALSE,
    do_double_rotation=TRUE,
    do_planar_fit=FALSE,
    do_flagging=TRUE, dir_blocked=c(0,0),
    do_SNDcorrection=TRUE,A=7/8,B=7/8,
    #do_WPLcorrection=FALSE,
    store=TRUE,format_out="txt",filename=NULL,
    meta=TRUE
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
    lint=time_averaging*60/(time_resolution) #length of averaging interval, i.e. number of measurements to be averaged
    #prepare output data
    cat("\n... allocate storage for output data ...")
    column_names=c("u_mean","v_mean","w_mean","Ts_mean","h2o_mean","co2_mean","ch4_mean",
                    "u_sd","v_sd","w_sd","Ts_sd","h2o_sd","co2_sd","ch4_sd",
                    "wd","ws",
                    "tke","ustar","L","zeta",
                    "cov_uw","cov_vw","cov_uv","cov_wTs","cov_vTs","cov_h2ow","cov_co2w","cov_ch4w",
                    "cov_wT_snd",
                    "sh","lh","co2_flux","methane_flux",
                    "flag_all","flag_stationarity","flag_distortion","flag_w","flag_itc",
                    "rotation_angle1","rotation_angle2")
    out=array(NA,dim=c(nint,length(column_names)))
    out=as.data.frame(out)
    colnames(out)=column_names
    #despiking
    if (do_despiking==TRUE) {
        cat("\n... do despiking ...")
        u=despiking(u,c(despike_u[1],despike_u[2]),despike_u[3],despike_u[4],despike_v[5])
        v=despiking(v,c(despike_v[1],despike_v[2]),despike_v[3],despike_v[4],despike_v[5])
        w=despiking(w,c(despike_w[1],despike_w[2]),despike_w[3],despike_w[4],despike_w[5])
        temp=despiking(temp,c(despike_temp[1],despike_temp[2]),despike_temp[3],despike_temp[4],despike_temp[5])
        if (do_h2o) h2o=despiking(h2o,c(despike_h2o[1],despike_h2o[2]),despike_h2o[3],despike_h2o[4],despike_h2o[5])
        if (do_co2) co2=despiking(co2,c(despike_co2[1],despike_co2[2]),despike_co2[3],despike_co2[4],despike_co2[5])
        if (do_ch4) ch4=despiking(ch4,c(despike_ch4[1],despike_ch4[2]),despike_ch4[3],despike_ch4[4],despike_ch4[5])
    }
    #rotation pre-check
    if (do_double_rotation==TRUE & do_planar_fit==TRUE) warning("You have chosen two rotation types, but only one can be applied. Apply double rotation now.")
    if (do_double_rotation==FALSE & do_planar_fit==FALSE) warning("You have chosen no rotation type, so no rotation is applied to the data.")
    cat("\n... start loop over data: do double rotation and stationarity flagging (if requested) ...")
    for (i in 1:nint) {
        #cat(paste0("\n\t #index: ",i,"\t progress: ",round(i/nint*100,2)," %"))
        i1=(i-1)*lint+1
        i2=(i*lint)
        iselect=seq(i1,i2)
        #wind (before rotation, assumes that the sonic is oriented towards north as indicated on the instrument)
        out$ws[i]=mean(calc_windSpeed2D(u[iselect],v[iselect]),na.rm=T)
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
        if (do_flagging == TRUE) {
            out$flag_stationarity=flag_stationarity(temp[iselect],w[iselect],nsub=as.integer(lint/4))
        }
    }
    #planar fit
    if (do_planar_fit) {
        wind_rotated=rotate_planar(u,v,w)
        u=wind_rotated$u
        v=wind_rotated$v
        w=wind_rotated$w
        out$rotation_angle1=wind_rotated$alpha
        out$rotation_angle2=wind_rotated$beta
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
    u_avg=averaging(u,time_resolution,time_averaging*60)
    out$u_mean=u_avg$mean[[1]]
    out$u_sd=u_avg$sd[[1]]
    v_avg=averaging(v,time_resolution,time_averaging*60)
    out$v_mean=v_avg$mean[[1]]
    out$v_sd=v_avg$sd[[1]]
    w_avg=averaging(w,time_resolution,time_averaging*60)
    out$w_mean=w_avg$mean[[1]]
    out$w_sd=w_avg$sd[[1]]
    Ts_avg=averaging(temp,time_resolution,time_averaging*60)
    out$Ts_mean=Ts_avg$mean[[1]]
    out$Ts_sd=Ts_avg$sd[[1]]
    if (do_h2o) {
        h2o_avg=averaging(h2o,time_resolution,time_averaging*60)
        out$h2o_mean=h2o_avg$mean[[1]]
        out$h2o_sd=h2o_avg$sd[[1]]
    }
    if (do_co2) {
        co2_avg=averaging(co2,time_resolution,time_averaging*60)
        out$co2_mean=co2_avg$mean[[1]]
        out$co2_sd=co2_avg$sd[[1]]
    }
    if (do_ch4) {
        ch4_avg=averaging(ch4,time_resolution,time_averaging*60)
        out$ch4_mean=ch4_avg$mean[[1]]
        out$ch4_sd=ch4_avg$sd[[1]]
    }
    #flux calculation
    cat("\n... do flux calculation ...")
    out$cov_uw=averaging(u*w,time_resolution,time_averaging*60)$mean[[1]]
    out$cov_uv=averaging(u*v,time_resolution,time_averaging*60)$mean[[1]]
    out$cov_vw=averaging(v*w,time_resolution,time_averaging*60)$mean[[1]]
    out$cov_wTs=averaging(w*temp,time_resolution,time_averaging*60)$mean[[1]]
    out$cov_vTs=averaging(v*temp,time_resolution,time_averaging*60)$mean[[1]]
    if (do_h2o) out$cov_h2ow=averaging(w*h2o,time_resolution,time_averaging*60)$mean[[1]]
    if (do_co2) out$cov_co2w=averaging(w*co2,time_resolution,time_averaging*60)$mean[[1]]
    if (do_ch4) out$cov_ch4w=averaging(w*ch4,time_resolution,time_averaging*60)$mean[[1]]
    #SND correction
    if (do_SNDcorrection==TRUE) {
        if (do_h2o == FALSE) {
            out$cov_wT_snd=SNDcorrection(out$Ts_mean,out$u_mean,out$v_mean,out$cov_uw,out$cov_vw,out$cov_wTs,NULL,A,B)
        } else {
            out$cov_wT_snd=SNDcorrection(out$Ts_mean,out$u_mean,out$v_mean,out$cov_uw,out$cov_vw,out$cov_wTs,out$cov_qw,A,B)
        }
        out$sh=cov2sh(out$cov_wT_snd)
    } else {
        out$sh=cov2sh(out$cov_wTs)
    }
    #WPL correction
    #if (do_WPLcorrection) {
    #    if (do_h2o) cov_h2ow_wpl=WPLcorrection(h2o,w,temp)
    #    if (do_co2) cov_co2w_wpl=WPLcorrection(co2,w,temp)
    #    if (do_ch4) cov_ch4w_wpl=WPLcorrection(ch4,w,temp)
    #}
    if (do_h2o) out$lh=cov2lh(out$cov_h2ow)
    #if (do_h2o & do_WPLcorrection) out$lh=cov2lh(cov_h2ow_wpl)
    if (do_co2) out$co2_flux=cov2cf(out$cov_co2w)
    #if (do_co2 & do_WPLcorrection) out$co2_flux=cov2cf(cov_co2w_wpl)
    #calculate turbulence statistics
    out$tke=calc_tke(out$u_sd,out$v_sd,out$w_sd)
    out$ustar=calc_ustar(out$cov_uw,out$cov_vw)
    out$L=calc_L(out$ustar,out$Ts_mean,out$cov_wTs)
    out$zeta=calc_zeta(measurement_height,out$L)
    #flagging: the other flags
    if (do_flagging==TRUE) {
        out$flag_itc=flag_most(out$w_sd,out$ustar,out$zeta)
        out$flag_w=flag_w(out$w_mean)
        out$flag_all=pmax(out$flag_stationarity,out$flag_distortion,out$flag_itc,out$flag_w,na.rm=T)
    }
    #------------------------------------------------
    #store post-processed data
    if (is.null(filename)) {
        systime=Sys.time()
        systime_string=format(systime,"%F_%H%M%S",tz="utc")
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
            "\ndo_planar_fit: ", do_planar_fit,
            "\ndo_flagging: ", do_flagging,
            "\ndo_SNDcorrection: ", do_SNDcorrection
        )
        cat(paste0("\n\n",meta,"\n\n"))
        if (store==TRUE) writeLines(meta,paste0("metadata_ec-processing_Reddy_",systime_string,".txt"))
    }
    return(out)
}
