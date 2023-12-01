#' Flux-Footprint Parametrization (FFP) according to Kljun et al., 2015
#'
#'@description Calculates Flux-Footprint Parametrization (FFP) according to Kljun et al., 2015
#'@param zm measurement height [m]
#'@param umean mean horizontal wind speed [m/s] (alternatively you can also use z0)
#'@param h boundary-layer height [m]
#'@param L Obukhov length [m]
#'@param sigmav standard deviation of crosswind [m/s]
#'@param ustar friction velocity [m/s]
#'@param z0 roughness length [m] (either umean or z0 have to be given)
#'@param nres resolution (default is nres=1000)
#'
#'@return 
#'@export
#'
#'@examples
#'ffp=calc_flux_footprint(zm=20,z0=0.01,h=200,L=-100,sigmav=0.6,ustar=0.4,contours=0.8)
#'
calc_flux_footprint = function(zm, umean=NA, h, L, sigmav, ustar, z0=NA,contours=seq(0.9,0.1,-0.1),nres=1000) {
    #fitting parameter for crosswind-integrated footprint, see (Kljun et al., 2015) eq. 17
    a=1.452
    b=-1.991
    c=1.462
    d=0.136
    xstarmax=-c/b+d #eq. 20
    #fitting parameter for crosswind footprint, see (Kljun et al., 2015) eq. 19
    ac=2.17
    bc=1.66
    cc=20
    #prepare non-dimensional upwind distance
    xstar=seq(d,30,length=nres)
    #calculate crosswind-integrated footprint
    fstar=a*(xstar-d)^b*exp(-c/(xstar-d)) #eq. 14
    #calculate standard deviation of crosswind distance
    sigmay_star=ac*sqrt(bc*xstar^2/(1+cc*xstar)) #eq. 18
    #calculate real scale footprint and maximum location
    if (!is.na(umean)) { #use umean if given
        aux1=zm/(1-zm/h)*umean/ustar*karman()
        x=xstar*aux1 #eq. 21 for all xstar
        fy_mean=fstar/aux1 #eq. 8 inverted
        xmax=xstarmax*aux1
    } else if (!is.na(z0)) {
        #uses eq. 5, stability correction function for diabatic winprofiles, following Hogstroem, 1996
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
        print("ERROR: You have to know either umean or z0.")
    }
    print(aux2)
    #calculate real scale sigmay
    ps1=min(1,abs(1/(zm/L))*10^-5 + ifelse(L<=0,0.8,0.55))
    sigmay=sigmay_star/ps1*zm*sigmav/ustar #eq. 13 inverted
    #calculate real scale f(x,y)
    dx=x[3]-x[2]
    print(length(x))
    print(dx)
    ypos=seq(0,length(x)/2*dx*1.5,dx)
    #print(ypos)
    fpos=matrix(NA,nrow=length(fy_mean),ncol=length(ypos))
    for (i in 1:length(fy_mean)) {
        fpos[i,]=fy_mean[i]/(sqrt(2*pi)*sigmay[i])*exp(-ypos^2/(2*sigmay[i]^2)) #eq. 10
    }
    #print(fy_mean)
    #footprint 2d
    n1=1
    n2=length(ypos)
    nf=nrow(fpos)
    print(n1)
    print(n2)
    fmat=array(NA,dim=c(nf,2*n2-1))
    fmat[,1:(n2-1)]=fpos[,n2:2]
    fmat[,n2:(2*n2-1)]=fpos
    #print(fmat)
    nx=length(x)
    y=c(-rev(ypos),ypos[2:n2])
    ny=length(y)
    xmat=matrix(rep(y,nx),nrow=ny,ncol=nx,byrow=T)
    ymat=matrix(rep(x,ny),nrow=ny,ncol=nx)
    #footprint contours
    nc=length(contours)
    fmat_tmp=rev(sort(c(fmat)))
    fmat_tmp=fmat_tmp[!is.na(fmat_tmp)]
    fint=cumsum(fmat_tmp)*dx^2
    ffp_cont=list()
    for (i in 1:nc) {
        cont=contourLines(x,y,fmat,levels=contours[i])
        print(str(cont))
        ffp_cont$xcont[[i]]=cont[[1]]$x
        ffp_cont$ycont[[i]]=cont[[1]]$y
    }
    #output
    ffp$xmax=xmax
    ffp$x=x
    ffp$y=y
    ffp$x2d=xmat
    ffp$y2d=ymat
    ffp$f2d=fmat
    ffp$xcontour=ffp_cont$xcont
    ffp$ycontour=ffp_cont$ycont
}