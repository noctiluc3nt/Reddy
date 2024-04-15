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
#'a=rnorm(100)
#'b=rnorm(100)
#'qa_ab=calc_quadrant_analysis(a,b)
#'
calc_quadrant_analysis=function(xval,yval,do_normalization=TRUE,hole_sizes=seq(0,10)) {
    covariance_total=cov(xval,yval,use="pairwise.complete.obs")
    correlation_total=cor(xval,yval,use="pairwise.complete.obs")
    if (do_normalization==TRUE) {
        xval = (xval-mean(xval,na.rm=T))/sd(xval,na.rm=T)
        yval = (yval-mean(yval,na.rm=T))/sd(yval,na.rm=T)
    }
    product_total=xval*yval
    nh = length(hole_sizes)
    occurrence = array(NA,dim=c(4,nh))
    product = array(NA,dim=c(4,nh))
    covariance = array(NA,dim=c(4,nh))
    for (i in 1:nh) {
        #q1
        is.q1 = (xval>(hole_sizes[i]/xval) & yval>(-hole_sizes[i]/yval))
        occurrence[1,i] = sum(is.q1,na.rm=TRUE)
        product[1,i] = mean(xval[is.q1] * yval[is.q1],na.rm=T)
        covariance[1,i] = cov(xval[is.q1],yval[is.q1])
        #q2
        is.q2 = (xval<(-hole_sizes[i]/xval) & yval>(-hole_sizes[i]/xval))
        occurrence[2,i] = sum(is.q2,na.rm=TRUE)
        product[2,i] = mean(xval[is.q2] * yval[is.q2],na.rm=T)
        covariance[2,i] = cov(xval[is.q2],yval[is.q2])
        #q3
        is.q3 = (xval<(hole_sizes[i]/xval) & yval<(hole_sizes[i]/yval))
        occurrence[3,i] = sum(is.q3,na.rm=TRUE)
        product[3,i] = mean(xval[is.q3] * yval[is.q3],na.rm=T)
        covariance[3,i] = cov(xval[is.q3],yval[is.q3])
        #q4
        is.q4 = (xval>(-hole_sizes[i]/xval) & yval<(-hole_sizes[i]/xval))
        occurrence[4,i] = sum(is.q4,na.rm=TRUE)
        product[4,i] = mean(xval[is.q4] * yval[is.q4],rm=T)
        covariance[4,i] = cov(xval[is.q4],yval[is.q4])
    }
    return(list("hole_sizes"=hole_sizes,
            "occurrence"=occurrence,
            "product"=product,
            "covariance"=covariance,
            "covariance_total"=covariance_total,
            "correlation_total"=correlation_total,
            "product_total"=product_total,
            "meta"="Output format: rows represent the quadrants Q1, Q2, Q3, Q4 -- columns represent selected hole sizes"))
}

#' Plotting Quadrant Analysis
#'
#'@description Calculates occurrence fraction and strength of the four quadrants
#'@param xval values of x variable (vector)
#'@param yval values of y variable (vector)
#'@param do_normalization should the values be normalized, i.e. (x-mean(x))/sd(x)? default: TRUE
#'@param hole_sizes vector containing desired hole sizes (integers >= 0), default: c(1,2)
#'@param contours vector containing levels of contour lines for 2d kernel densoty estimation, default: contours=10^(-3:3)
#'@param print_fit should the fit summary from the linear regression be printed? default: TRUE
#'@param ... arguments passed to plot function
#'@return 
#'@export
#'
#'@importFrom MASS kde2d
#' 
#'@examples
#'a=rnorm(100)
#'b=rnorm(100)
#'plot_quadrant_analysis(a,b)
#'
plot_quadrant_analysis=function(xval,yval,do_normalization=TRUE,hole_sizes=c(1,2),contours=10^(-3:3),print_fit=TRUE,...) {
    if (do_normalization) {
        xval=(xval-mean(xval,na.rm=T))/sd(xval)
        yval=(yval-mean(yval,na.rm=T))/sd(yval)
    }
    if (!exists("pch")) pch = 20
    if (typeof(col)=="closure") col = rgb(0.6,0.6,0.6,0.1)
    plot(xval,yval,col=col,pch=pch,...)
    abline(h=0,lty=2,col=1,lwd=2)
    abline(v=0,lty=2,col=1,lwd=2)
    #linear regression
    fit=lm(yval~xval)
    if (print_fit==TRUE) print(summary(fit))
    abline(fit,lwd=3,col="red4")
    #2d kde
    nc=length(contours)
    lab=colorRampPalette(c("blue3","red3"), space = "Lab")
    kde=MASS::kde2d(xval,yval)
    contour(kde$x,kde$y,kde$z,levels=contours,col=lab(nc)[1:nc],add=TRUE,lwd=2,drawlabels=FALSE)
    #draw holes
    xs=seq(-10,10,0.1)
    for (h in hole_sizes) {
        points(xs,h/xs,col=1,type="l",lty=2,lwd=1,pch=pch)
        points(xs,-h/xs,col=1,type="l",lty=2,lwd=1,pch=pch)
    }
}


