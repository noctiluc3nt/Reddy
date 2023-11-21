#' Invariant analysis of Reynolds stress tensor
#'
#'@description Invariant analysis of Reynolds stress tensor, calculation of anisotropy
#'@param a11 R11 element of Reynolds stress tensor (scalar or vector)
#'@param a12 R12 element of Reynolds stress tensor (scalar or vector)
#'@param a13 R13 element of Reynolds stress tensor (scalar or vector)
#'@param a22 R22 element of Reynolds stress tensor (scalar or vector)
#'@param a23 R23 element of Reynolds stress tensor (scalar or vector)
#'@param a33 R33 element of Reynolds stress tensor (scalar or vector)
#'@return list containing eta, xi, xb, yb (eta,xi are the coordinates of the Lumley triangle and xb,yb the coordinates of the barycentric map)
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
	#symmetry
	a21=a12
	a31=a13
	a32=a23
	for (i in 1:n) {
		rey=matrix(c(a11[i]^2,a12[i],a13[i],a21[i], a22[i]^2,a23[i],a31[i],a32[i],a33[i]^2),nrow=3,ncol=3)
		#cat("\n",i,rey)
		if (!any(is.na(rey))) {	
			if (sum(diag(rey))!=0) {	
				B = rey/sum(diag(rey))-1/3*delta
				#diagonalize the anisotropy matrix (calculate the eigenvalues)			
				ev=eigen(B)$values
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
			}
		}
	}
	return(list("eta"=eta,"xi"=xi,"xb"=xb,"yb"=yb))
}