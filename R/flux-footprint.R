#' Flux-Footprint Parametrization (FFP) according to Kljun et al., 2015
#'
#'@description Calculates the Flux-Footprint Parametrization (FFP) according to Kljun et al., 2015
#'@param zm measurement height [m]
#'@param ws_mean mean horizontal wind speed [m/s] (alternatively you can also use \code{z0})
#'@param wd_mean mean wind direction [deg] (used to rotate flux footprint, optional)
#'@param L Obukhov length [m]
#'@param v_sd standard deviation of crosswind [m/s]
#'@param ustar friction velocity [m/s]
#'@param z0 roughness length [m] (either \code{ws_mean} or \code{z0} have to be given)
#'@param blh boundary-layer height [m]
#'@param contours which contour lines should be calculated? default: \code{contours=seq(0.9,0.1,-0.1)}
#'@param nres resolution (default: \code{nres=1000})
#'@param plot logical, should the flux footprint be plotted? default \code{plot=TRUE}
#'
#'@return list containing all relevant flux footprint information
#'@export
#'
#'@importFrom omnibus rotateMatrix
#'
#'@examples
#'#unrotated (i.e. if wind direction not given):
#'ffp=calc_flux_footprint_Kljun2015(zm=20,ws_mean=2,blh=200,L=-1.5,
#'      v_sd=0.6,ustar=0.4,contours=0.8)
#'#rotated (i.e. given wind direction):
#'ffp=calc_flux_footprint_Kljun2015(zm=20,ws_mean=2,wd_mean=80,blh=200,L=-1.5,
#'      v_sd=0.6,ustar=0.4,contours=0.8)
#'
calc_flux_footprint_Kljun2015 = function(zm,ws_mean=NA,wd_mean=NA,L,v_sd,ustar,z0=NA,blh,contours=seq(0.9,0.1,-0.1),nres=1000,plot=TRUE) {
    #fitting parameters for crosswind-integrated footprint, see (Kljun et al., 2015) eq. 17
    a=1.452
    b=-1.991
    c=1.462
    d=0.136
    xstarmax=-c/b+d #eq. 20
    #fitting parameters for crosswind footprint, see (Kljun et al., 2015) eq. 19
    ac=2.17
    bc=1.66
    cc=20
    #prepare non-dimensional upwind distance
    xstar=seq(d,30,length=nres)[-1]
    #calculate crosswind-integrated footprint
    fstar=a*(xstar-d)^b*exp(-c/(xstar-d)) #eq. 14
    #calculate standard deviation of crosswind distance
    sigmay_star=ac*sqrt(bc*xstar^2/(1+cc*xstar)) #eq. 18
    #calculate real scale footprint and maximum location
    if (!is.na(ws_mean)) { #use ws_mean if given
        aux1=zm/(1-zm/blh)*ws_mean/ustar*karman()
        x=xstar*aux1 #eq. 21 for all xstar
        fy_mean=fstar/aux1 #eq. 8 inverted
        xmax=xstarmax*aux1
    } else if (!is.na(z0)) {
        #uses eq. 5, stability correction function for diabatic wind profiles following Hogstroem, 1996
        if (L<=0) { #unstable or neutral conditions
            xi=(1-19*zm/L)^0.25
            psim=log((1+xi^2)/2) + 2*log((1+xi)/2) - 2*atan(xi)+pi/2 
        } else { #stable conditions
            psim=-5.3*zm/L
        }
        aux2=zm/(1-zm/blh)*(log(zm/z0)-psim)
        x=xstar*aux2 #eq. 22
        xmax=xstarmax*aux2
        if ((log(zm/z0)-psim)>0) {
            fy_mean=fstar/aux2 #eq. 9
        } else {
            error = -1
        }
    } else {
        stop("You have to know either ws_mean or z0.")
    }
    #calculate real scale y_sd
    ps1=min(1,abs(1/(zm/L))*1E-5 + ifelse(L<=0,0.8,0.55))
    sigmay=sigmay_star/ps1*zm*v_sd/ustar #eq. 13 inverted
    #calculate real scale f(x,y)
    dx=x[3]-x[2]
    ypos=x #seq(0,length(x)/2*dx*1.5,dx)
    fpos=matrix(NA,nrow=length(fy_mean),ncol=length(ypos))
    for (i in 1:length(fy_mean)) {
        fpos[i,]=fy_mean[i]/(sqrt(2*pi)*sigmay[i])*exp(-ypos^2/(2*sigmay[i]^2)) #eq. 10
    }
    #footprint 2d
    n1=1
    n2=length(ypos)
    nf=nrow(fpos)
    fmat=array(NA,dim=c(nf,2*n2-1))
    fmat[,1:(n2-1)]=fpos[,n2:2] #symmetry
    fmat[,n2:(2*n2-1)]=fpos
    nx=length(x)
    y=c(-rev(ypos),ypos[2:n2])
    ny=length(y)
    #xmat=matrix(rep(x,ny),nrow=ny,ncol=nx,byrow=TRUE)
    #ymat=matrix(rep(y,nx),nrow=ny,ncol=nx)
    #cut coordinates/fmat
    ny_mid=round(ny/2)
    yselect=(ny_mid-nres/2+1):(ny_mid+nres/2)
    fmat=fmat[1:(nres/2),yselect]
    fmat=rbind(fmat*0,fmat)
    x=seq(min(-x),max(x),length.out=dim(fmat)[1])
    y=x
    #rotate flux footprint
    if (!is.na(wd_mean)) {
        #rotate area (2d data) using polar coordinates
        #angle=atan2(ymat,xmat)
        #dist=sqrt(xmat^2+ymat^2)
        #xmat_rot=dist*sin(wd-angle)
        #ymat_rot=dist*cos(wd-angle)
        #rotate ffp directly using omnibus package
        wd_mean=ifelse(wd_mean<180,wd_mean+180,wd_mean-180)
        fmat=omnibus::rotateMatrix(fmat,wd_mean+90)
        fmat[which(is.na(fmat),arr.ind=TRUE)]=0
        #rotate contours
        #for (i in 1:length(contours)) {  
        #    wd=wd_mean*pi/180
        #    angle=atan2(ffp_cont$ycont[[i]],ffp_cont$xcont[[i]])
        #    dist=sqrt(ffp_cont$xcont[[i]]^2+ffp_cont$ycont[[i]]^2)
        #    xcont_rot[[i]]=dist*sin(wd-angle)
        #    ycont_rot[[i]]=dist*cos(wd-angle)
        #}
    }
    #footprint contours
    fcont=get_contours_from_f2d(x,y,fmat,contours=contours)
    xcont=fcont$xcont
    ycont=fcont$ycont
    #nc=length(contours)
    #fsort=rev(sort(c(fmat)))
    #fsort=fsort[!is.na(fsort)]
    #fint=cumsum(fsort)*dx^2
    #xcont=list()
    #ycont=list()
    #for (i in 1:nc) {
    #    fdiff=abs(fint-contours[i])
    #    ind=which.min(fdiff)
    #    fr=fsort[ind]
    #    cont=contourLines(x,y,fmat,levels=fr)
    #    xcont[[i]]=cont[[1]]$x
    #    ycont[[i]]=cont[[1]]$y
    #}
    #output
    ffp=list()
    ffp$xmax=xmax
    ffp$x=x
    ffp$y=y
    ffp$fy_mean=fy_mean
    ffp$f2d=fmat
    ffp$xcontour=xcont
    ffp$ycontour=ycont
    ffp$contour_levels=contours
    tryCatch({
        if (plot==TRUE) plot_flux_footprint(ffp)
    }, warning = function(w) {
        message("An error occurred when plotting the flux footprint.")
    })
    return(ffp)
}


#' Flux-Footprint Calculation according to Korman and Meixner, 2001
#'
#'@description Calculates the Flux-Footprint Parametrization (FFP) according to Korman and Meixner, 2001
#'@param zm measurement height [m]
#'@param ws_mean mean horizontal wind speed [m/s] (alternatively you can also use \code{z0})
#'@param wd_mean mean wind direction [deg] (used to rotate flux footprint, optional)
#'@param L Obukhov length [m]
#'@param v_sd standard deviation of crosswind [m/s]
#'@param ustar friction velocity [m/s]
#'@param z0 roughness length [m] (either \code{ws_mean} or \code{z0} have to be given)
#'@param contours which contour lines should be calculated? default: \code{contours=seq(0.9,0.1,-0.1)}
#'@param nres domain size (default: \code{nres=1000} to get a domain ranging from -500 to 500)
#'@param dx resolution (default: \code{dx=1})
#'@param plot logical, should the flux footprint be plotted? default \code{plot=TRUE}
#'
#'@return list containing all relevant flux footprint information
#'@export
#'
#'@examples
#'#unrotated (i.e. if wind direction not given):
#'ffp=calc_flux_footprint_KM2001(zm=20,ws_mean=2,L=-1.5,
#'      v_sd=0.6,ustar=0.4,z0=0.1,contours=0.8)
#'#rotated (i.e. given wind direction):
#'ffp=calc_flux_footprint_KM2001(zm=20,ws_mean=2,wd_mean=80,L=-1.5,
#'      v_sd=0.6,ustar=0.4,z0=0.1,contours=0.8)
#'
calc_flux_footprint_KM2001 = function(zm,ws_mean=NA,wd_mean=NA,L,v_sd,ustar,z0,contours=seq(0.9,0.1,-0.1),nres=1000,dx=1,plot=TRUE) {
    #ctodo heck all same length
    #init + create domain
    ngrid=as.integer(nres/2)
    x=seq(-ngrid,ngrid,dx)
    y=seq(-ngrid,ngrid,dx)
    nx=length(x)
    grid=expand.grid(x,y)
    x2d=matrix(grid[,1],nrow=nx,byrow=T)
    y2d=matrix(grid[,2],nrow=nx,byrow=T)
    f2d=array(NA,dim=c(nx,nx))
    zeta=zm/L
    if (zeta>0) { #stable
        phi=1+5*zeta #eq 33
        if (!is.na(ws_mean)) { #real wind
            u=ws_mean
        } else { #MOST windprofile
            u=ustar/karman()*(log(zm/z0)+5*zeta) #eq 31 with 35 for psi
        }
        m=ustar/karman()*phi/u #eq 36
        n=1/phi #eq 36
    } else if (zeta==0) { #neutral
        phi=1
        if (!is.na(ws_mean)) { #real wind
            u=ws_mean
        } else { #MOST windprofile
            u=ustar/karman()*log(zm/z0)
        }
        m=ustar/karman()/u
        n=1
    } else { #unstable
        zeta2=(1-16*zeta)^0.25 #eq 35 below
        psi=-2*log((1+zeta2)/2)-log((1+zeta2^2)/2)+2*atan(zeta2)-pi/2 #eq 35
        phi=(1-16*zeta)^(-0.25) #eq 33
        if (!is.na(ws_mean)) { #real wind
            u=ws_mean
        } else { #MOST windprofile
            u=ustar/karman()*(log(zm/z0)+psi) #eq 31
        }
        m=ustar/karman()*phi/u
        n=(1-24*zeta)/(1-16*zeta) #eq 36
    }
    r=2+m-n #table 1 / shape factor
    mu=(1+m)/r #table 1 continued
    U=u/(zm^m) #eq 11 / A1
    K=karman()*ustar*zm/phi/zm^n #A2
    xi=U*zm^r/(r^2*K) #eq 19 / flux length scale
    gamma_mu=ifelse(mu==0,0,gamma(mu))
    gamma_1r=ifelse(1/r==0,0,gamma(1/r))      
    #parameters from ART footprint tool from Neftel et al, 2007
    A=1+mu #A9
    B=U*(zm^r)/(r^2*K) #A10
    C=(B^mu)/gamma_mu #A11
    D=v_sd*gamma_1r*(r^2*K/U)^(m/r)/(gamma_mu*U) #A12
    E=(r-m)/r #A13
    #wind direction
    wd=ifelse(!is.na(wd_mean),wd_mean,0) #(wd_mean-90)%%360 #ifelse(wd<180,wd+180,wd-180)
    x_tmp=x2d*cos(wd*pi/180)+y2d*sin(wd*pi/180)
    y_tmp=-x2d*sin(wd*pi/180)+y2d*cos(wd*pi/180)
    x_mask=(x_tmp>0)
    x_masked=x_mask*x_tmp
    #calculation        
    h0=D*x_masked^E
    f2d=1/(sqrt(2*pi)*h0)*exp(-y_tmp^2/(2*h0^2))*C*x_masked^(-A)*exp(-B/x_masked) #eq 2 in Neftel et al, 2007
    f2d[which(is.na(f2d))]=0
    #output
    ffp=list()
    ffp$x=x
    ffp$y=y
    ffp$f2d=f2d
    tryCatch({
        #get contours
        fcont=get_contours_from_f2d(x,y,f2d,contours=contours)
        xcont=fcont$xcont
        ycont=fcont$ycont
        ffp$xcontour=xcont
        ffp$ycontour=ycont
        ffp$contour_levels=contours
    }, error = function(w) {
        message("An error occurred when calculating the contours.")
    })
    tryCatch({
        if (plot==TRUE) plot_flux_footprint(ffp)
    }, warning = function(w) {
        message("An error occurred when plotting the flux footprint.")
    })
    return(ffp)
}


#' Get contours from 2D flux footprint matrix
#'
#'@description Calculates contours for given flux footprint
#'@param x x-coordinate (vector)
#'@param y y-coordinate (vector)
#'@param fmat 2D flux footprint (matrix/array with the dimensions matching nx x ny) 
#'@param contours which contours? default \code{contours=seq(0,0.9,0.1)}
#' 
#'@return list with x- and y-coordinates of the contours
#'@export
#' 
get_contours_from_f2d=function(x,y,fmat,contours=seq(0,0.9,0.1)) {
    #footprint contours
    nc=length(contours)
    fsort=rev(sort(c(fmat)))
    fsort=fsort[!is.na(fsort)]
    dx=x[3]-x[2]
    fint=cumsum(fsort)*dx^2
    xcont=list()
    ycont=list()
    for (i in 1:nc) {
        fdiff=abs(fint-contours[i])
        ind=which.min(fdiff)
        fr=fsort[ind]
        cont=contourLines(x,y,fmat,levels=fr)
        xcont[[i]]=cont[[1]]$x
        ycont[[i]]=cont[[1]]$y
    }
    return(list(xcont,ycont))
}

#' Plot Flux-Footprint
#'
#'@description Plots Flux-Footprint Parametrization (FFP)
#'@param ffp an object returned from \code{calc_flux_footprint_[method]}
#'@param levels levels used for filled contour plot of footprint, default \code{levels=c(0,10^seq(-6,-3,0.1))}
#'@param mode can be either \code{mode="distance"} for plotting footprint relative to station location in cartesian coordinates or \code{mode="lonlat"} for plotting in (lon,lat)-ccordinates
#'@param ... paraemters passed to image.plot function 
#'
#'@return no return
#'@importFrom grDevices colorRampPalette contourLines rgb
#'@importFrom graphics .filled.contour abline arrows axis
#'@importFrom stats approx ccf cor cov lm median quantile sd
#'@importFrom fields image.plot
#'@export
#'
#'@examples
#'ffp=calc_flux_footprint_Kljun2015(zm=5,ws_mean=5,blh=700,L=-1.3,v_sd=1.2,ustar=0.35)
#'plot_flux_footprint(ffp)
#' 
plot_flux_footprint = function(ffp,levels=c(0,10^seq(-6,-3,0.1)),mode="distance",...) {
    #filled contour plot with contour lines
    #lab=colorRampPalette(c("white","blue3","yellow","orange","red3"), space = "Lab")
    #nlev=length(levels)
    #plot(NA,xlim=xlim,ylim=ylim,main="2D Flux Footprint",xlab="x [m]",ylab="y [m]")
    #.filled.contour(ffp$x2d[1,],ffp$y2d[,1],ffp$f2d,levels=levels,col=lab(nlev))
    if (mode=="distance") {
        if (!exists("xlim")) xlim=c(-200,200)
        if (!exists("ylim")) ylim=c(-200,200)
        #plot crosswind-integrated footprint
        #tryCatch({
        #    plot(ffp$x,ffp$fy_mean,type="l",xlim=xlim,lwd=2,xlab="x [m]",ylab="crosswind-integrated footprint",main="Crosswind-Integrated Flux Footprint")
        #    abline(v=ffp$xmax,col=2,lwd=2,lty=2)
        #    legend("topright",legend="footprint peak location",col=2,lwd=2,lty=2)
        #})
        #fields::image.plot(ffp$x2d,ffp$y2d,ffp$f2d*100,xlim=xlim,ylim=ylim,xlab="x [m]",ylab="y [m]",main="Flux Footprint")
        collab=colorRampPalette(c("white","lightblue","blue3","red3"), space = "Lab")
        fields::image.plot(ffp$x,ffp$y,ffp$f2d*100,xlim=xlim,ylim=ylim,xlab="x [m]",ylab="y [m]",main="Flux Footprint",col=collab(64))
        points(0,0,pch=20)
        tryCatch({
            for (i in 1:length(ffp$xcontour)) {
                lines(ffp$xcontour[[i]],ffp$ycontour[[i]],type="l",lwd=1)
            }
        })
    } else if (mode=="lonlat") {
        if (is.null(ffp$x2d_earth)) {
            message("The flux footprint has so far been calculated only in cartesian coordinates --
            you need to use locate_flux_footprint to transform it to (lon,lat)-coordinates")
        } else {
            if (!exists("xlim")) xlim=range(ffp$x2d_earth)
            if (!exists("ylim")) ylim=range(ffp$y2d_earth)
            fields::image.plot(ffp$x2d_earth,ffp$y2d_earth,ffp$f2d*100,xlab="lon",ylab="lat",main="Flux Footprint",xlim=xlim,ylim=ylim)
            tryCatch({
                for (i in 1:length(ffp$xcontour)) {
                    lines(ffp$xcontour_earth[[i]],ffp$ycontour_earth[[i]],type="l",lwd=1)
                }
            })
        }
    }
    #3d perspective plot
    #nmid=as.integer(length(ffp$x2d[,1])/2)
    #xselect=(nmid-200):(nmid+200)
    #persp(ffp$x2d[1,xselect],ffp$y2d[xselect,1],ffp$f2d[xselect,xselect],main="2d flux footprint as 3d plot",xlab=xlab,ylab=ylab,zlab="footprint",theta = 30, phi = 30, expand = 0.5, col = "lightblue",ltheta = 120, shade = 0.2,nticks=5)
}



#' Transform flux footprint from (x,y)-coordinates to (lon,lat)-coordinates through given station location
#'
#'@description 
#'@param ffp ffp object returned by \code{calc_flux_footprint_[method]}
#'@param lon_station lon of station location
#'@param lat_station lat of station location
#'
#'@return ffp object which contains ffp also in lon-lat coordinates (with extension _earth)
#'@export
#'
#'@examples
#'lon1=7.527061462
#'lat1=60.59384155
#'ffp=calc_flux_footprint_Kljun2015(zm=20,ws_mean=2,blh=200,L=-1.5,v_sd=0.6,ustar=0.4,contours=0.8)
#'ffp=locate_flux_footprint(ffp,lon1,lat1)
#'
locate_flux_footprint = function(ffp,lon_station,lat_station) {
    x2lon=function(xmat,lon_s=lon_station,lat_s=lat_station) {
        return(lon_s + xmat/R_earth()*180/pi*cos(lat_s*pi/180))
    }
    y2lat=function(ymat,lat_s=lat_station) {
        return(lat_s + ymat/R_earth()*180/pi)
    }
    ffp$x_earth = x2lon(ffp$x)
    ffp$y_earth = y2lat(ffp$y)
    ffp$xcontour_earth = lapply(ffp$xcontour,x2lon)
    ffp$ycontour_earth = lapply(ffp$ycontour,y2lat)
    return(ffp)
}

#' Calculate Flux Footprint Climatology
#'
#'@description Calculates the Flux-Footprint Parametrization (FFP) according to Kljun et al., 2015
#'@param zm measurement height [m]
#'@param ws_mean mean horizontal wind speed [m/s] (alternatively you can also use \code{z0})
#'@param wd_mean mean wind direction [deg] (used to rotate flux footprint, optional)
#'@param blh boundary-layer height [m]
#'@param L Obukhov length [m]
#'@param v_sd standard deviation of crosswind [m/s]
#'@param ustar friction velocity [m/s]
#'@param z0 roughness length [m] (either \code{ws_mean} or \code{z0} have to be given)
#'@param contours which contour lines should be calculated? default: \code{contours=seq(0.9,0.1,-0.1)}
#'@param nres resolution (scalar) (default: \code{nres=1000})
#'@param method method used to calculate FFP: can be either \code{method="Kljun2015"} (default) or \code{method="KM2001"}
#'@param plot logical, should the flux footprint be plotted? default \code{plot=TRUE}
#'@param ... paraemters passed to image.plot function
#'
#'@return list containing all relevant flux footprint information
#'@export
#'
#'@examples
#'nit=2
#'#ffp=calc_flux_footprint_climatology(zm=20,ws_mean=rep(2,nit),blh=rep(200,nit),L=rep(-1.5,nit),
#'#          v_sd=rep(0.6,nit),ustar=rep(0.4,nit),contours=0.8,wd_mean=c(50,240)) #todo
#'
calc_flux_footprint_climatology = function(zm,ws_mean=NA,wd_mean=NA,L,v_sd,ustar,z0=NA,blh=NA,contours=seq(0.9,0.1,-0.1),nres=1000,method="Kljun2015",plot=TRUE,...) {
    n=length(ws_mean)
    ffp_clim=array(0,dim=c(nres+1,nres+1))
    ncount=0
    for (i in 1:n) {
        cat("progress: ",round(i/n*100))
        tryCatch({
            if (method=="Kljun2015") {
                ffp_tmp=calc_flux_footprint_Kljun2015(zm=zm,ws_mean=ws_mean[i],wd_mean=wd_mean[i],blh=blh[i],L=L[i],v_sd=v_sd[i],ustar=ustar[i],z0=z0[i],nres=nres,plot=FALSE)
            } else if (method=="KM2001") {
                ffp_tmp=calc_flux_footprint_KM2001(zm=zm,ws_mean=ws_mean[i],wd_mean=wd_mean[i],L=L[i],v_sd=v_sd[i],ustar=ustar[i],z0=z0[i],nres=nres,plot=FALSE)
            } else {
                warning("The requested flux footprint calculation method is not implemented: method has to be either Kljun2015 or KM2001.")
            }
            f2d_tmp=ffp_tmp$f2d
            f2d_tmp[which(is.na(f2d_tmp),arr.ind=TRUE)]=0
            ffp_clim=ffp_clim+f2d_tmp
            ncount=ncount+1
        },
        error = function(e) {
            message("an error occurred.")
        })
    }
    x=ffp_tmp$x
    y=ffp_tmp$y
    fmat=ffp_clim/ncount
    #output
    ffp=list()
    #ffp$xmax=xmax
    ffp$x=x
    ffp$y=y
    tryCatch({
        #get contours
        fcont=get_contours_from_f2d(x,y,fmat,contours=contours)
        xcont=fcont$xcont
        ycont=fcont$ycont
        ffp$xcontour=xcont
        ffp$ycontour=ycont
        ffp$contour_levels=contours
    }, error = function(w) {
        message("An error occurred when calculating the contours.")
    })
    tryCatch({
        if (plot==TRUE) plot_flux_footprint(ffp,...)
    }, warning = function(w) {
        message("An error occurred when plotting the flux footprint.")
    })
    return(ffp)
}
