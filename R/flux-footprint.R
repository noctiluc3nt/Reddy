#' Flux-Footprint Parametrization (FFP) according to Kljun et al., 2015
#'
#'@description Calculates the Flux-Footprint Parametrization (FFP) according to Kljun et al., 2015
#'@param zm measurement height [m]
#'@param u_mean mean horizontal wind speed [m/s] (alternatively you can also use \code{z0})
#'@param h boundary-layer height [m]
#'@param L Obukhov length [m]
#'@param v_sd standard deviation of crosswind [m/s]
#'@param ustar friction velocity [m/s]
#'@param z0 roughness length [m] (either \code{u_mean} or \code{z0} have to be given)
#'@param contours which contour lines should be calculated? default: \code{contours=seq(0.9,0.1,-0.1)}
#'@param nres resolution (default: \code{nres=1000})
#'@param plot logical, should the flux footprint be plotted? default \code{plot=TRUE}
#'
#'@return list containing all relevant flux footprint information
#'@export
#'
#'@examples
#'ffp=calc_flux_footprint(zm=20,u_mean=2,h=200,L=-1.5,v_sd=0.6,ustar=0.4,contours=0.8)
#'
calc_flux_footprint = function(zm, u_mean=NA, h, L, v_sd, ustar, z0=NA,contours=seq(0.9,0.1,-0.1),nres=1000,plot=TRUE) {
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
    if (!is.na(u_mean)) { #use u_mean if given
        aux1=zm/(1-zm/h)*u_mean/ustar*karman()
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
        aux2=zm/(1-zm/h)*(log(zm/z0)-psim)
        x=xstar*aux2 #eq. 22
        xmax=xstarmax*aux2
        if ((log(zm/z0)-psim)>0) {
            fy_mean=fstar/aux2 #eq. 9
        } else {
            error = -1
        }
    } else {
        stop("You have to know either u_mean or z0.")
    }
    #calculate real scale sigmay
    ps1=min(1,abs(1/(zm/L))*1E-5 + ifelse(L<=0,0.8,0.55))
    sigmay=sigmay_star/ps1*zm*v_sd/ustar #eq. 13 inverted
    #calculate real scale f(x,y)
    dx=x[3]-x[2]
    ypos=seq(0,length(x)/2*dx*1.5,dx)
    fpos=matrix(NA,nrow=length(fy_mean),ncol=length(ypos))
    for (i in 1:length(fy_mean)) {
        fpos[i,]=fy_mean[i]/(sqrt(2*pi)*sigmay[i])*exp(-ypos^2/(2*sigmay[i]^2)) #eq. 10
    }
    #footprint 2d
    n1=1
    n2=length(ypos)
    nf=nrow(fpos)
    fmat=array(NA,dim=c(nf,2*n2-1))
    fmat[,1:(n2-1)]=fpos[,n2:2]
    fmat[,n2:(2*n2-1)]=fpos
    nx=length(x)
    y=c(-rev(ypos),ypos[2:n2])
    ny=length(y)
    xmat=matrix(rep(x,ny),nrow=ny,ncol=nx,byrow=TRUE)
    ymat=matrix(rep(y,nx),nrow=ny,ncol=nx)
    #footprint contours
    nc=length(contours)
    fsort=rev(sort(c(fmat)))
    fsort=fsort[!is.na(fsort)]
    fint=cumsum(fsort)*dx^2
    ffp_cont=list()
    for (i in 1:nc) {
        fdiff=abs(fint-contours[i])
        ind=which.min(fdiff)
        fr=fsort[ind]
        cont=contourLines(x,y,fmat,levels=fr)
        ffp_cont$xcont[[i]]=cont[[1]]$x
        ffp_cont$ycont[[i]]=cont[[1]]$y
    }
    #output
    ffp=list()
    ffp$xmax=xmax
    ffp$x=x
    ffp$fy_mean=fy_mean
    ffp$x2d=xmat
    ffp$y2d=ymat
    ffp$f2d=fmat
    ffp$xcontour=ffp_cont$xcont
    ffp$ycontour=ffp_cont$ycont
    ffp$contour_levels=contours
    if (plot==TRUE) plot_flux_footprint(ffp)
    return(ffp)
}


#' Plot Flux-Footprint
#'
#'@description Plots Flux-Footprint Parametrization (FFP) according to Kljun et al., 2015
#'@param ffp an object returned from \code{calc_flux_footprint}
#'@param levels levels used for filled contour plot of footprint, default \code{levels=c(0,10^seq(-6,-3,0.1))}
#' 
#'@return no return
#'@importFrom grDevices colorRampPalette contourLines rgb
#'@importFrom graphics .filled.contour abline arrows axis
#'@importFrom stats approx ccf cor cov lm median quantile sd
#'@export
#'
#'@examples
#'ffp=calc_flux_footprint(zm=5,u_mean=5,h=700,L=-1.3,v_sd=1.2,ustar=0.35)
#'plot_flux_footprint(ffp)
#' 
plot_flux_footprint = function(ffp,levels=c(0,10^seq(-6,-3,0.1))) {
    if (!exists("xlim")) xlim=c(0,400)
    if (!exists("ylim")) ylim=c(-250,250)
    #plot crosswind-integrated footprint
    plot(ffp$x,ffp$fy_mean,type="l",xlim=xlim,lwd=2,xlab="x [m]",ylab="crosswind-integrated footprint",main="Crosswind-Integrated Flux Footprint")
    abline(v=ffp$xmax,col=2,lwd=2,lty=2)
    legend("topright",legend="footprint peak location",col=2,lwd=2,lty=2)
    #filled contour plot with contour lines
    lab=colorRampPalette(c("white","blue3","yellow","orange","red3"), space = "Lab")
    nlev=length(levels)
    plot(NA,xlim=xlim,ylim=ylim,main="2D Flux Footprint",xlab="x [m]",ylab="y [m]")
    .filled.contour(ffp$x2d[1,],ffp$y2d[,1],ffp$f2d,levels=levels,col=lab(nlev))
    for (i in 1:length(ffp$xcontour)) {
        lines(ffp$xcontour[[i]],ffp$ycontour[[i]],type="l",lwd=1)
    }
    #3d perspective plot
    #nmid=as.integer(length(ffp$x2d[,1])/2)
    #xselect=(nmid-200):(nmid+200)
    #persp(ffp$x2d[1,xselect],ffp$y2d[xselect,1],ffp$f2d[xselect,xselect],main="2d flux footprint as 3d plot",xlab=xlab,ylab=ylab,zlab="footprint",theta = 30, phi = 30, expand = 0.5, col = "lightblue",ltheta = 120, shade = 0.2,nticks=5)
}