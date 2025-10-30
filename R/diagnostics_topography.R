
#' Calculates surface roughness length z0 from friction velocity using the simple estimate from Charnock, 1955
#'
#'@description Calculates surface roughness z0 from friction velocity using the simple estimate from Charnock, 1955: z0 = alpha*ustar^2/g with alpha=0.016 and g=9.81 m/s^2
#'@param ustar friction velocity [m/s]
#'
#'@return surface roughness length [m]
#'@export
#'
#'@examples
#'ustar2z0(0.2)
#'
ustar2z0 = function(ustar) {
	return(alpha()*ustar^2/g())
}


#' Calculates the slope-based direction-dependent terrain parameter Sx that indicates the degree of shelter or exposure of a point of interest from a given digital elevation model, according to Winstral and Marks, 2002
#'
#'@description Calculates the slope-based direction-dependent terrain parameter Sx from a digital elevation model
#'@param xi x-coodinate of point of interest
#'@param yi y-coodinate of point of interest
#'@param elevi elevation of point of interest
#'@param x x-coordinates (2-dimensional field)
#'@param y y-coordinates (2-dimensional field)
#'@param elev elevations (2-dimensional field)
#'@param wd wind direction
#'@param dmax maximum fetch distance considered
#'@param mode either \code{mode="cartesian"} (default) or \code{mode="lonlat"} depending on whether (x, y) are in cartesian coordinates or representing longitude, latitude 
#'
#'@return Sx (with Sx>0 indicating shelter and Sx<0 exposure)
#'@export
#'
#'
calc_Sx = function(xi,yi,elevi,x,y,elev,wd,dmax=1000,mode="cartesian") {
	if (mode=="cartesian") {
        angles=atan2((y-yi),(x-xi))
        dists=sqrt((x-xi)^2+(y-yi)^2)
    } else if (mode=="lonlat") {
        dists=calc_distance(x,y,xi,yi)
        #todo angles
    } else {
        warning("The mode has to be either cartesian or lonlat.")
    }
    cond=(dists<dmax & angles==wd)
    Sx=max(atan2(elev[cond]-elevi,dists[cond]))
    return(Sx)
}


#' Calculates the topographic position index TPI
#'
#'@description Calculates the topographic position index TPI (direction-independent)
#'@param xi x-coodinate of point of interest
#'@param yi y-coodinate of point of interest
#'@param elevi elevation of point of interest
#'@param x x-coordinates (2-dimensional field)
#'@param y y-coordinates (2-dimensional field)
#'@param elev elevations (2-dimensional field)
#'@param dmax maximum fetch distance considered
#'@param mode either \code{mode="cartesian"} (default) or \code{mode="lonlat"} depending on whether (x, y) are in cartesian coordinates or representing longitude, latitude 
#'
#'@return TPI (with TPI>0 indicating exposure and TPI<0 shelter, i.e. elevation lower than everage)
#'@export
#'
#'
calc_TPI = function(xi,yi,elevi,x,y,elev,dmax=1000,mode="cartesian") {
	if (mode=="cartesian") {
        dists=sqrt((x-xi)^2+(y-yi)^2)
    } else if (mode=="lonlat") {
        dists=calc_distance(x,y,xi,yi)
    } else {
        warning("The mode has to be either cartesian or lonlat.")
    }
    cond=(dists<dmax)
    TPI=elevi-mean(elev[cond],na.rm=T)
    return(TPI)
}


#' Calculates the 'deviation from mean elevation' (DEV)
#'
#'@description Calculates the 'deviation from mean elevation' (DEV) (direction-independent) as TPI/SD (normalized)
#'@param xi x-coodinate of point of interest
#'@param yi y-coodinate of point of interest
#'@param elevi elevation of point of interest
#'@param x x-coordinates (2-dimensional field)
#'@param y y-coordinates (2-dimensional field)
#'@param elev elevations (2-dimensional field)
#'@param dmax maximum fetch distance considered
#'@param mode either \code{mode="cartesian"} (default) or \code{mode="lonlat"} depending on whether (x, y) are in cartesian coordinates or representing longitude, latitude 
#'
#'@return DEV (with DEV>0 indicating exposure and DEV<0 shelter, i.e. elevation lower than everage)
#'@export
#'
#'
calc_DEV = function(xi,yi,elevi,x,y,elev,dmax=1000,mode="cartesian") {
    TPI=calc_TPI(xi,yi,elevi,x,y,elev,dmax=1000,mode="cartesian")
    if (mode=="cartesian") {
        dists=sqrt((x-xi)^2+(y-yi)^2)
    } else if (mode=="lonlat") {
        dists=calc_distance(x,y,xi,yi)
    } else {
        warning("The mode has to be either cartesian or lonlat.")
    }
    cond=(dists<dmax)
    SD=sd(elev[cond],na.rm=T)
    DEV=TPI/SD
    return(DEV)
}
