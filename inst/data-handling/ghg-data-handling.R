#' unzip directory containing .ghg files
#'
#'@description unzips .ghg files from input directory (\code{dir_in}) to output directory (\code{dir_out})
#'@param dir_in input directory path (character)
#'@param dir_out output directory path (character)
#'
#'@return no return
#'@export
#'
unzip_dir = function(dir_in,dir_out) {
	files=list.files(dir_in,pattern="ghg",full.names=TRUE)
    nf=length(files)
    for (i in 1:nf) unzip(files[i],exdir=dir_out)
    #maybe add: remove .status and .metadata files in dir_out, such that only .data files remain
}


#' read ghg data (unzipped)
#'
#'@description reads unzipped ghg-files (.data) and transforms them into a dataframe containing: time, u, v, w, Ts, q, co2, pressure
#'@param file file (as character, extension .data) after unzipping (with \code{unzip_ghg()})
#'@param do_h2o logical, read H2O as well? default \code{do_h2o=FALSE}
#'@param do_co2 logical, read CO2 as well? default \code{do_co2=FALSE}
#'@param do_ch4 logical, read CH4 as well? default \code{do_ch4=FALSE}
#'
#'@return dataframe, with: time [beginn of measurement], Ts [K], u [m/s], v [m/s], w [m/s], pressure [Pa], h2o -- if measured [original unit, usually milli mol], co2 -- if measured [original unit, usually mikro mol], ch4 -- if measured [original unit, usually mikro mol]
#'@export
#'
read_ghg_data = function(file,do_h2o=FALSE,do_co2=FALSE,do_ch4=FALSE) {
	dat=read.table(file,sep="\t",skip=7,header=T)
    tryCatch({
	    date=as.character(dat$Date)
	    time=substr(as.character(dat$Time),1,12)
        out=data.frame("Ts"=dat$SOS..m.s., #sos
			"u"=dat$U..m.s.,
			"v"=dat$V..m.s.,
			"w"=dat$W..m.s.,
            "pressure"=dat$Total.Pressure..kPa.*1000)
        if (do_h2o) out$h2o=dat$H2O.dry.mmol.mol. #milli mol (name depends on device)
	    if (do_co2) out$co2=dat$CO2.dry.umol.mol. #mikro mol		
   	    if (do_ch4) out$ch4=dat$CH4..umol.mol. #mikro mol		
	    out$Ts=sos2Ts(out$Ts) #speed of sound to sonic temperature
        out$date=date
        out$time=time
        return(out)
    }, error = function(e) {
        message("An error occurred while reading the file -- this could e.g. be caused by an incomplete file.")
    })
}
