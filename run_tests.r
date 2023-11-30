source("./R/anisotropy.R")
source("./R/mrd.R")


### example MRD ###
series=c(1,3,2,5,1,2,1,3)
calc_mrd(series)

set.seed(5)
series=rnorm(2^10)
mrd_test=calc_mrd(series)
plot_mrd(mrd_test)

### example anisotropy ###
calc_anisotropy(1,0,0,1,0,1) #isotropic (i.e., yb=0.8660254 corresponding to sqrt(3)/2)
calc_anisotropy(1,0,0.5,1,0,1) #some anisotropy