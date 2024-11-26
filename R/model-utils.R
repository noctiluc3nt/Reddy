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
		field_out=dat
   		for (i in 2:n) {
        	out[,,i]=dat[,,i]-dat[,,i-1] #deaccumulate hourly
    	}
	}
    return(out*factor) #time unit conversion
}


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
