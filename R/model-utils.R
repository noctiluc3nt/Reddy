#' deaccumulation 
#'
#'@description hourly deaccumulation, e.g. for fluxes from model output
#'@param dat vector (with dimension time) or array (with dimension x, y, time)
#'@param factor factor for unit and sign conversion, default: \code{factor = -1/3600} for converting hour to seconds and adapting the sign convention
#'@return vector or array hourly deaccumulated (same dimension as input)
#'@export
#'
deaccumulate1h=function(dat,factor=-1/3600) {
	if (is.vector(dat)) {
		n=length(dat)
		out=dat
    	for (i in 2:n) {
        	out[i]=dat[i]-dat[i-1] #deaccumulate hourly
    	}
	} else {
		dims=dim(dat)
		n=dims[3]
		out=dat
   		for (i in 2:n) {
        	out[,,i]=dat[,,i]-dat[,,i-1] #deaccumulate hourly
    	}
	}
    return(out*factor) #time unit conversion
}

################ coordinates #################################
#' Converts hybrid (terrain-following) sigma levels to physical heights 
#'
#'@description Converts hybrid (terrain-following) sigma levels to physical heights
#'@param hybrid scalar or vector, hybrid sigma levels
#'@param Tv virtual temperature
#'@return hybrid levels converted to physical height [m]
#'@export
#'
#'@examples
#'sigma2height(0.1)
#'sigma2height(0.1,288)
#'
sigma2height=function(hybrid,Tv=273.15) {
    return(-Rd()*Tv/g()*log(hybrid))
}

#' Converts height ot pressure  
#'
#'@description Converts physical height to pressure based on given surface pressure and temperature
#'@param height height/elevation [m]
#'@param pres0 surface pressure [Pa], default \code{pres0=101315}
#'@param temp0 reference temperature [k], default \code{temp0=273.15}
#'@return pressure [Pa] (at the given height levels)
#'@export
#'
#'@examples
#'height2pres(c(10,100,1000),pres=99900)
#'
height2pres=function(height,pres0=101315,temp0=273.15) {
    return(pres0*(1-height*0.0065/temp0)^5.255)
}

###################### derivatives ###########################

#' df_dx
#'
#'@description Calculates x-derivative for equidistant grid
#'@param fld input field with dimension (x,y)
#'@param xres resolution in x-direction
#'@return x-derivative of fld (same dimensions)
#'@export
#'
#'@examples
#'set.seed(5)
#'field=matrix(rnorm(16),ncol=4)
#'df_dx(field,10)
#'
df_dx=function(fld,xres=1) {
    nx=dim(fld)[1]
    ny=dim(fld)[2]
    df_dx=array(NA,dim=c(nx,ny)) 
    df_dx[1,]=(fld[2,]-fld[1,])/xres #left boundary
    df_dx[nx,]=(fld[nx,]-fld[nx-1,])/xres #right boundary
    for (i in 2:(nx-1)) {
        df_dx[i,]=(fld[i+1,]-fld[i-1,])/(2*xres) #interior
    }
    return(df_dx)
}

#' df_dy
#'
#'@description Calculates y-derivative for equidistant grid
#'@param fld input field with dimension (x,y)
#'@param yres resolution in y-direction
#'@return y-derivative of fld (same dimensions)
#'@export
#'
#'@examples
#'set.seed(5)
#'field=matrix(rnorm(16),ncol=4)
#'df_dy(field,10)
#'
df_dy=function(fld,yres=1) {
    nx=dim(fld)[1]
    ny=dim(fld)[2]
    df_dy=array(NA,dim=c(nx,ny)) 
    df_dy[,1]=(fld[,2]-fld[,1])/yres #lower boundary
    df_dy[,ny]=(fld[,ny]-fld[,ny-1])/yres #upper boundary
    for (i in 2:(ny-1)) {
        df_dy[,i]=(fld[,i+1]-fld[,i-1])/(2*yres) #interior
    }
    return(df_dy)
}

#################### distance ######################################
#' Calculates distance between two points on a sphere given in (lon, lat)-coordinates
#'
#'@description Calculates distance between two (lon, lat)-points on a spheroid
#'@param lon1 longitude location 1 [deg]
#'@param lat1 latitude location 1 [deg]
#'@param lon2 longitude location 2 [deg]
#'@param lat2 latitude location 2 [deg]
#'@return distance between two points [m]
#'@export
#'
#'@examples
#'calc_distance(8,60,8,61)
#'
calc_distance=function(lon1,lat1,lon2,lat2) {
    lon1=lon1*pi/180
    lat1=lat1*pi/180
    lon2=lon2*pi/180
    lat2=lat2*pi/180
    dlon=lon2-lon1
    dlat=lat2-lat1
    h1=sin(dlat/2)^2+cos(lat1)*cos(lat2)*sin(dlon/2)^2
    h2=2*atan2(sqrt(h1),sqrt(1-h1))
    return(R_earth()*h2)
}

#' Find closest grid point
#'
#'@description Finds the closest grid point from a given point (lon_loc, lat_loc) to a given grid (lons,lats) 
#'@param lons longitudes of a grid (vector or matrix) [deg]
#'@param lats latitudes of a grid (vector or matrix) [deg]
#'@param lon_loc longitude of the location [deg]
#'@param lat_loc latitude of the location [deg]
#'@return array index of lons, lats that represents the closest grid point to lon_loc, lat_loc
#'@export
#'
find_closest_grid_point = function(lons,lats,lon_loc,lat_loc) {
    dist=sqrt((lons-lon_loc)^2 + (lats-lat_loc)^2)
    ind=which(dist==min(dist),arr.ind=T)
    return(ind)
}
