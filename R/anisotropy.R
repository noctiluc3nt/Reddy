#' Invariant analysis of Reynolds stress tensor
#'
#'@description Invariant analysis of Reynolds stress tensor, calculation of Lumley and barycentric map coordinates and anisotropy
#'@param a11 R11 element of Reynolds stress tensor: u_sd^2   (scalar or vector)
#'@param a12 R12 element of Reynolds stress tensor: cov(u,v) (scalar or vector)
#'@param a13 R13 element of Reynolds stress tensor: cov(u,w) (scalar or vector)
#'@param a22 R22 element of Reynolds stress tensor: v_sd^2   (scalar or vector)
#'@param a23 R23 element of Reynolds stress tensor: cov(v,w) (scalar or vector)
#'@param a33 R33 element of Reynolds stress tensor: w_sd^2   (scalar or vector)
#'@return list containing xb, yb, eta, xi, all eigenvalues and eigenvectors (eta, xi are the coordinates of the Lumley triangle and xb, yb the coordinates of the barycentric map)
#'@export
#'
#'@examples
#'calc_anisotropy(1,0,0,1,0,1) #isotropic
#'calc_anisotropy(1,0,1,1,0,1) #some anisotropy
#'
calc_anisotropy = function(a11,a12,a13,a22,a23,a33) {
	n=length(a11)
	#unit matrix / Kronecker delta
	delta=matrix(c(1,0,0,0,1,0,0,0,1), ncol=3,nrow=3)
	#initialize
	eta=numeric(n)
	xi=numeric(n)
	xb=numeric(n)
	yb=numeric(n)
	evals=array(NA,dim=c(n,3))
	evecs=array(NA,dim=c(n,3,3))
	#symmetry
	a21=a12
	a31=a13
	a32=a23
	for (i in 1:n) {
		rey=matrix(c(a11[i],a12[i],a13[i],a21[i], a22[i],a23[i],a31[i],a32[i],a33[i]),nrow=3,ncol=3)
		#cat("\n",i,rey)
		if (!any(is.na(rey))) {	
			if (sum(diag(rey))!=0) {	
				B = rey/sum(diag(rey))-1/3*delta
				#diagonalize the anisotropy matrix (calculate the eigenvalues)
				inv=eigen(B) #invariant analysis		
				ev=inv$values #eigenvalues
				Bdiag=ev*delta
				evs_sort=rev(sort(ev))
				eta[i] = (1/3*(evs_sort[1]^2 + evs_sort[1]*evs_sort[2] + evs_sort[2]^2))^(1/2)
				dummyX=evs_sort[1]*evs_sort[2]*(evs_sort[1]+evs_sort[2])
				xi[i] = -sign(dummyX)*(1/2*abs(dummyX))^(1/3)
				#baricentric map
				C1c = evs_sort[1]-evs_sort[2]
				C2c = 2*(evs_sort[2]-evs_sort[3])
				C3c = 3*(evs_sort[3]) + 1
				#coordinates
				xb[i] = C1c + 0.5*C3c
				yb[i] = C3c*sqrt(3)/2
				#store
				evals[i,]=ev
				evecs[i,,]=inv$vectors
			}
		}
	}
	return(list("xb"=xb,"yb"=yb,"eta"=eta,"xi"=xi,"eigenvalues"=evals,"eigenvectors"=evecs))
}


#' Plot in barycentric map
#'
#'@description Plots (xb, yb) from invariant analysis of Reynolds stress tensor (calc_anisotropy) in barycentric map
#'@param xb xb coordinate (e.g., from calc_anisotropy)
#'@param yb yb coordinate (e.g., from calc_anisotropy)
#'@param contours vector containing levels of contour lines for 2d kernel densoty estimation, default: contours=c(5,10,20)
#'
#'@return plots (xb, yb) in barycentric map with 2d kernel density estimation (no return)
#'@export
#' 
#'@importFrom MASS kde2d
#'
#'@examples
#'nm=100
#'example1=calc_anisotropy(rep(1,nm),rep(0,nm),runif(nm,0,1),rep(1,nm),rep(0,nm),runif(nm,1,1.5))
#'plot_barycentric_map(example1$xb,example1$yb)
#'
plot_barycentric_map = function(xb,yb,contours=c(5,10,20)) {
	if (!exists("pch")) pch = 20
    if (typeof(col)=="closure") col = rgb(0,0,0,0.4)
    plot(xb,yb,pch=pch,col=col,xlim=c(0,1),ylim=c(0,sqrt(3)/2),asp=1,xlab="x",ylab="x",main="Barycentric Map")
    segments(0,0,1,0,lwd=2)
	segments(0,0,0.5,sqrt(3)/2,lwd=2)
	segments(1,0,0.5,sqrt(3)/2,lwd=2)
	points(0.5,sqrt(3)/6,pch=3,lwd=2)
    #2d kde
    nc=length(contours)
    lab=colorRampPalette(c("blue3","red3"), space = "Lab")
    kde=MASS::kde2d(xb,yb)
    contour(kde$x,kde$y,kde$z,levels=contours,col=lab(nc)[1:nc],add=TRUE,lwd=2,drawlabels=FALSE)
}