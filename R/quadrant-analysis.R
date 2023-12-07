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
#' a=rnorm(100)
#'
calc_quadrant_analysis=function(xval,yval,do_normalization=TRUE,hole_sizes=seq(0,10)) {
    if (do_normalization==TRUE) {
        xval = (xval-mean(xval,na.rm=T))/sd(xval,na.rm=T)
        yval = (yval-mean(yval,na.rm=T))/sd(yval,na.rm=T)
    }
    nh = length(hole_sizes)
    occurrence = array(NA,dim=c(4,nh))
    product = array(NA,dim=c(4,nh))
    covariance = array(NA,dim=c(4,nh))
    for (i in 1:nh) {
        #q1
        is.q1 = (xval>hole_sizes[i]/yval & yval>hole_sizes[i]/xval)
        occurrence[1,i] = sum(is.q1,na.rm=TRUE)
        product[1,i] = mean(xval[is.q1] * yval[is.q1],na.rm=T)
        covariance[1,i] = cov(xval[is.q1],yval[is.q1])
        #q2
        is.q2 = (xval<(-hole_sizes[i]/yval) & yval>(-hole_sizes[i]/xval))
        occurrence[2,i] = sum(is.q2,na.rm=TRUE)
        product[2,i] = mean(xval[is.q2] * yval[is.q2],na.rm=T)
        covariance[2,i] = cov(xval[is.q2],yval[is.q2])
        #q3
        is.q3 = (xval<hole_sizes[i]/yval & yval<hole_sizes[i]/xval)
        occurrence[3,i] = sum(is.q3,na.rm=TRUE)
        product[3,i] = mean(xval[is.q3] * yval[is.q3],na.rm=T)
        covariance[3,i] = cov(xval[is.q3],yval[is.q3])
        #q4
        is.q4 = (xval>(-hole_sizes[i]/yval) & yval<(-hole_sizes[i]/xval))
        occurrence[4,i] = sum(is.q4,na.rm=TRUE)
        product[4,i] = mean(xval[is.q4] * yval[is.q4],rm=T)
        covariance[4,i] = cov(xval[is.q4],yval[is.q4])
    }
    return(list(hole_sizes,occurrence,product,covariance))
}

#' Plotting Quadrant Analysis TODO
#'
#'@description Calculates Occurrence Fraction and Strength of the four Quadrants
#'@param xval values of x variable (vector)
#'@param yval values of y variable (vector)
#'@param do_normalization should the values be normalized? i.e. (x-mean(x))/sd(x)
#'@param hole_sizes vector containing desired hole sizes (integers >= 0)
#'@param ... arguments passed to plot function
#'@return 
#'@export
#'
#'@examples
#' #TODO
#'
calc_quadrant_analysis=function(xval,yval,do_normalization=TRUE,hole_sizes=seq(0,10),...) {
    if (do_normalization) {
        xval=(xval-mean(xval,na.rm=T))/sd(xval)
        yval=(yval-mean(yval,na.rm=T))/sd(yval)
    }
    #TODO plot with contour lines and hole lines
    plot(xval,yval,...)
}


