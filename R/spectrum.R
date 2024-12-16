#' Spectrum of timeseries by wrapping rbase::spectrum()
#'
#'@description Calculates and plots the averaged turbulence spectrum (as wrapper of rbase::spectrum)
#'@param ts timeseries
#'@param nbins number of bins used to average the spectrum, default \code{nbins=100}
#'@param plot should the spectrum be plotted? default \code{plot=TRUE}
#'
#'@return binned spectrum
#'@export
#'
#'@examples
#'set.seed(5)
#'ts=rnorm(1000)
#'calc_spectrum(ts,nbins=100,plot=FALSE)
#'
calc_spectrum = function(ts,nbins=100,plot=TRUE) {
	s=spectrum(ts,plot=FALSE)
    bins=seq(log(min(s$freq,na.rm=TRUE)),log(max(s$freq,na.rm=TRUE)),length.out=nbins)
    sbin=binning(s$spec,s$freq,10^bins)
    if (plot==TRUE) {
        plot(bins[2:nbins],log(sbin[,2]),pch=20,main="Spectrum",xlab="frequency",ylab="spectrum",xlim=range(bins[2:nbins][!is.na(sbin[,2])],na.rm=TRUE),xaxt="n",yaxt="n")
        axis(1,at=(-10):0,10^(-10:0))
        axis(2,at=(-10):10,10^(-10:10))
        fit=lm(log(sbin[,2]) ~ bins[2:nbins])
        print(summary(fit))
        #par(new=TRUE)
        #plot(bins,bins*(-5/3),type="l",col=4,lwd=2,lty=2,xlab="",ylab="",xaxt="n",yaxt="n",ylim=range(log(sbin[,2]),na.rm=TRUE))
        #legend("bottomleft",legend="-5/3 slope",col=4,lwd=2,lty=2)
    }
    return(cbind(bins[2:nbins],sbin[,2]))
}


#' Frequency spectrum (1D)
#'
#'@description Calculates and plots turbulence spectrum (in time) calculated using FFT (and optionally bins it)
#'@param ts timeseries of the variable for which the spectrum should be calculated
#'@param tres time resolution [s] of the given timeseries, default \code{tres=0.05} for 20 Hz
#'@param nbins number of bins used to average the spectrum, default \code{nbins=NULL}, i.e. no further binning is applied (means number of bins equals half of length of input timeseries)
#'@param method method used to calculate the spectrum, can be either FFT (fast Fourier transform) or DCT (discrete cosine transform), default FFT
#'@param na.rm should NA values be removed from the timeseries? default \code{na.rm=TRUE}
#'@param plot should the spectrum be plotted? default \code{plot=TRUE}
#'@param ... further arguments passed to plot function
#'
#'@importFrom gsignal dct
#'
#'@return binned frequency spectrum from 1D FFT
#'@export
#'
#'@examples
#'set.seed(5)
#'ts=rnorm(1000)
#'calc_spectrum1D(ts) #no binning
#'calc_spectrum1D(ts,nbins=100) #binning
#'
calc_spectrum1D = function(ts,tres=0.05,nbins=NULL,method="fft",na.rm=TRUE,plot=TRUE,...) {
    if (na.rm==TRUE) ts=ts[!is.na(ts)]
    nt=length(ts) #length of time series
    #sr=1/tres #sampling rate
    #t=seq(0,(nt-1)*tres,tres) #time vector
    if (method=="fft" | method=="FFT") {  #fft and frequencies
        spec=fft(ts)
        ssb=spec[1:(nt/2)] #single side band
        ssb[2:(nt/2)]=2*ssb[2:(nt/2)]
        spec=abs(ssb)*2/nt
        freq=seq(0,(nt/2)-1)*tres/nt #frequencies
        out=data.frame("frequency"=freq,"spectrum"=spec)
    } else if (method=="dct" | method=="DCT") { #dct and frequencies
        spec=abs(dct(ts))
        freq=seq(0,(nt)-1)*tres/nt #frequencies
        out=data.frame("frequency"=freq,"spectrum"=spec)
    }
    #binning
    #if (nbins>nt) warning("You request more bins than there are measurements available.")
    if (!is.null(nbins)) {  
        fbins=exp(seq(log(min(freq[-1])),log(max(freq)),length.out=nbins))
        xbinned=binning(spec,freq,fbins)
        fmid=(fbins[2:nbins]+fbins[1:(nbins-1)])/2
        out=data.frame("frequency"=fmid,"spectrum"=xbinned[,2])
    }
    if (plot==TRUE) {
        plot(out$frequency,out$spectrum,pch=20,log="xy",xlab="frequency [1/s]",...)
        #points(out$frequency,out$frequency^(-5/3),type="l",lty=2)
        #fit=lm(log(out$spectrum) ~ log(out$frequency))
        #print(summary(fit))
        #abline(exp(fit$coefficients[1]),exp(fit$coefficients[2]),col=2,lwd=2)
        #points(freq,freq^fit$coefficients[2]+exp(fit$coefficients[1]),col=2,type="l")
    }
    return(out)
}


#' Spatial spectrum (2D)
#'
#'@description Calculates and plots turbulence spectrum (in space) calculated using FFT or DCT (and optionally bins it)
#'@param field two-dimensional input field
#'@param xres spatial resolution in x-direction
#'@param yres spatial resolution in y-direction, default \code{yres=NULL} meaning that the field is equidistant and \code{yres=xres} is used
#'@param method method used to calculate the spectrum, can be either FFT (fast Fourier transform) or DCT (discrete cosine transform), default FFT
#'@param nbins number of bins used to average the spectrum, default \code{nbins=NULL}, i.e. no further binning is applied (means number of bins equals half of length of input timeseries)
#'@param plot should the spectrum be plotted? default \code{plot=TRUE}
#'@param ... further arguments passed to plot function
#'
#'@importFrom gsignal dct2 fftshift
#'
#'@return binned wavenumber spectrum from 2D FFT or DCT
#'@export
#'
#'@examples
#'set.seed(5)
#'field=matrix(rnorm(10000),nrow=100)
#'calc_spectrum2D(field,xres=100) #equidistant grid, no binning, fft
#'calc_spectrum2D(field,xres=100,yres=200) #non-equidistant grid, no binning, fft
#'calc_spectrum2D(field,xres=100,nbins=1000) #equidistant grid, binning, fft
#'calc_spectrum2D(field,xres=100,nbins=1000,method="dct") #equidistant grid, binning, dct
#'
calc_spectrum2D = function(field,xres=1000,yres=NULL,nbins=NULL,method="fft",plot=TRUE,...) {
    nx=dim(field)[1]
    ny=dim(field)[2]
    if (is.null(yres)) yres=xres
    #calc 2D spectrum
    if (method=="fft" | method=="FFT") {
        kx=seq(-nx/2,nx/2-1)/(nx*xres) #fftfreq(nx,xres)
        ky=seq(-ny/2,ny/2-1)/(ny*yres) #fftfreq(ny,yres)
        k=calc_k2d(kx,ky)
        spec=abs(fft(field)) #fft() automatically calculates 2D fft for 2D input array
        spec=gsignal::fftshift(spec,MARGIN=c(1,2)) #fftshift
    } else if (method=="dct" | method=="DCT") {
        kx=seq(0,nx-1)/(nx*xres)
        ky=seq(0,ny-1)/(ny*yres)
        k=calc_k2d(kx,ky)
        spec=abs(gsignal::dct2(field))
    }
    #binning
    if (is.null(nbins)) nbins=min(nx,ny)
    kbins=exp(seq(log(min(k[which(k>0,arr.ind=TRUE)])),log(max(k[which(is.finite(k),arr.ind=TRUE)])),length.out=nbins))
    xbinned=binning(c(spec),c(k),kbins)
    kmid=(kbins[2:nbins]+kbins[1:(nbins-1)])/2
    out=data.frame("wavenumber"=kmid,"spectrum"=xbinned[,2]) 
    if (plot==TRUE) {
        plot(out$wavenumber,out$spectrum,pch=16,log="xy",xlab="wavenumber [1/m]",ylab="spectrum",...)
    }
    return(out)
}


#' FFT frequency
#'
#'@description Returns FFT sampling frequencies (R version of the function fft.fftfreq in numpy)
#'@param n length
#'@param res spatial resolution (default: \code{res=1})
#'
#'@return vector of length n containing FFT frequencies
#'@export
#'
#'@examples
#'fftfreq(10)
#'
fftfreq=function(n,res=1) {
    if (n%%2==0) {
        freq=c(seq(0,(n/2-1)),seq((-n/2),-1))/(res*n)
    } else {
        freq=c(seq(0,((n-1)/2)),seq(-(n-1)/2,-1))/(res*n)
    }
    return(freq)
}

#' Calculates 2d (horizontal) wavenumber matrix from kx, ky
#'
#'@description Calculates 2d (horizontal) wavenumber matrix from kx, ky
#'@param kx wavenumber in x-direction
#'@param ky wavenumber in y-direction
#'
#'@return total spatial (horizontal) wavenumber k
#'@export
#'
#'@examples
#'kx=c(1:10)/10
#'ky=c(1:8)/8
#'k=calc_k2d(kx,ky)
#'
calc_k2d=function(kx,ky) {
    nx=length(kx)
    ny=length(ky)
    mat=array(NA,dim=c(nx,ny))
    for (i in 1:nx) {
        for (j in 1:ny) {
            mat[i,j]=(kx[i])^2+(ky[j])^2
        }
    }
    return(sqrt(mat))
}


#' Helmholtz-Hodge decomposition
#'
#'@description Calculates Helmholtz-Hodge decomposition of horizontal wind using a spectral FFT-based method: decomposition of horizontal wind in rotational and divergent part: (u,v) = (u_rot,v_rot) + (u_div,v_div)
#'@details The implementation is based on the Python version from https://github.com/shixun22/helmholtz.
#'@param u zonal wind field with dimension (x,y)
#'@param v meridional wind field with dimension (x,y)
#'@param res spatial resolution (assuming equidistant grid)
#'
#'@return list containing u_div, v_div, u_rot, v_rot
#'@export
#'
#'@examples
#'set.seed(5)
#'u=matrix(rnorm(100),ncol=10)
#'v=matrix(rnorm(100),ncol=10)
#'hd=calc_helmholtz_decomposition(u,v,100)
#'
calc_helmholtz_decomposition = function(u,v,res=1) {
    nx=dim(u)[1]
    ny=dim(u)[2]
    #fft
    ufft=fft(u)
    vfft=fft(v)
    #wave numbers
    kx=fftfreq(nx,res)
    ky=fftfreq(ny,res)
    k2=calc_k2d(kx,ky)^2
    k2[1,1]=1
    kx2=matrix(rep(kx,ny),nrow=nx,ncol=ny,byrow=F)
    ky2=matrix(rep(ky,nx),nrow=nx,ncol=ny,byrow=T)
    #decomposition
    div=ufft*kx2 + vfft*ky2
    u_div=gsignal::ifft(div/k2*kx2)
    v_div=gsignal::ifft(div/k2*ky2)
    u_rot=u-u_div
    v_rot=v-v_div
    return(list("u_div"=u_div,"v_div"=v_div,"u_rot"=u_rot,"v_rot"=v_rot))
}