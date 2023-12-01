#' Calculating Coherent Structures following Quadrant Analysis
#'
#'@description Calculates Occurrence Fraction and Strength of the four Quadrants
#'@param xval values of x variable (vector)
#'@param yval values of y variable (vector)
#'@param do_normalization should the values be normalized? i.e. (x-mean(x))/sd(x)
#'@param hole_sizes vector containing desired hole sizes (integers >= 0)
#'@return list containing occurrence fraction and strength (calculated based on product and covariance) for all four quadrants (mathematical orientation)
#'@export
#'
#'@examples
#' #TODO
#'
calc_quadrant_analysis=function(xval,yval,do_normalization=TRUE,hole_sizes=seq(0,10)) {
    if (do_normalization) {
        xval=(xval-mean(xval,na.rm=T))/sd(xval)
        yval=(yval-mean(yval,na.rm=T))/sd(yval)
    }
    nh=length(hole_sizes)
    occurrence=array(NA,dim=c(4,nh))
    product=array(NA,dim=c(4,nh))
    cov=array(NA,dim=c(4,nh))
    for (i in 1:nh) {
        #q1
        is.q1=(xval>hole_sizes[i]/yval & yval>hole_sizes[i]/xval)
        occurrence[1,i] = sum(is.q1,na.rm=TRUE)
        product[1,i] = xval[is.q1] * yval[is.q2]
        cov[1,i] = cov(xval[is.q1],yval[is.q2])
        #q2
        is.q2=(xval<(-hole_sizes[i]/yval) & yval>(-hole_sizes[i]/xval))
        occurrence[1,i] = sum(is.q1,na.rm=TRUE)
        product[1,i] = xval[is.q1] * yval[is.q2]
        cov[1,i] = cov(xval[is.q1],yval[is.q2])
        #q3
        is.q3=(xval<hole_sizes[i]/yval & yval<hole_sizes[i]/xval)
        occurrence[1,i] = sum(is.q1,na.rm=TRUE)
        product[1,i] = xval[is.q1] * yval[is.q2]
        cov[1,i] = cov(xval[is.q1],yval[is.q2])
        #q4
        is.q4=(xval>(-hole_sizes[i]/yval) & yval<(-hole_sizes[i]/xval))
        occurrence[1,i] = sum(is.q1,na.rm=TRUE)
        product[1,i] = xval[is.q1] * yval[is.q2]
        cov[1,i] = cov(xval[is.q1],yval[is.q2])
    }
    out=list(
        "occurrence_fraction"=occurrence,
        "strength_product"=product,
        "strength_cov"=cov
    )
    return(out)
}

#' Plotting Quadrant Analysis TODO