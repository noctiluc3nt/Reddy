#' Calculating Coherent Structures following Quadrant Analysis
#'
#'@description Calculates occurrence fraction and strength of the four quadrants in the framework of quadrant analysis
#'@param xval values of x variable (vector)
#'@param yval values of y variable (vector)
#'@param do_normalization should the values be normalized? i.e. \code{(x-mean(x))/sd(x)}, default: \code{do_normalization=TRUE}
#'@param hole_sizes vector containing desired hole sizes (integers >= 0)
#'@param orient only relevant for exuberance and organization ratio: if down-gradient flux corresponds to positive values, use \code{orient="+"} (for sensible and latent heat flux), if down-gradient flux corresponds to negative values, use \code{orient="-"} (for momentum flux and CO2 flux)
#'@param plot logical, should the quadrant analysis be plotted? default \code{plot=TRUE}
#'@param ... arguments passed to \code{plot_quadrant_analysis}
#'
#'@return list containing occurrence fraction and strength (calculated based on product and covariance) for all four quadrants (mathematical orientation) as well as the therefrom derived measures exuberance and organization ratio, i.e. the ratio of the strength (or occurrence frequency, respectively) of disorganized to organized structures
#'@export
#'
#'@examples
#'a=rnorm(100)
#'b=rnorm(100)
#'qa_ab=calc_quadrant_analysis(a,b)
#'
calc_quadrant_analysis=function(xval,yval,do_normalization=TRUE,hole_sizes=seq(0,10),orient="+",plot=TRUE,...) {
    covariance_total=cov(xval,yval,use="pairwise.complete.obs")
    correlation_total=cor(xval,yval,use="pairwise.complete.obs")
    if (do_normalization==TRUE) {
        xval = (xval-mean(xval,na.rm=TRUE))/sd(xval,na.rm=TRUE)
        yval = (yval-mean(yval,na.rm=TRUE))/sd(yval,na.rm=TRUE)
    }
    nh = length(hole_sizes)
    occurrence = array(NA,dim=c(4,nh))
    product = array(NA,dim=c(4,nh))
    covariance = array(NA,dim=c(4,nh))
    for (i in 1:nh) {
        #q1
        is.q1 = (xval>(hole_sizes[i]/xval) & yval>(-hole_sizes[i]/yval))
        occurrence[1,i] = sum(is.q1,na.rm=TRUE)
        product[1,i] = mean(xval[is.q1] * yval[is.q1],na.rm=TRUE)
        covariance[1,i] = cov(xval[is.q1],yval[is.q1])
        #q2
        is.q2 = (xval<(-hole_sizes[i]/xval) & yval>(-hole_sizes[i]/xval))
        occurrence[2,i] = sum(is.q2,na.rm=TRUE)
        product[2,i] = mean(xval[is.q2] * yval[is.q2],na.rm=TRUE)
        covariance[2,i] = cov(xval[is.q2],yval[is.q2])
        #q3
        is.q3 = (xval<(hole_sizes[i]/xval) & yval<(hole_sizes[i]/yval))
        occurrence[3,i] = sum(is.q3,na.rm=TRUE)
        product[3,i] = mean(xval[is.q3] * yval[is.q3],na.rm=TRUE)
        covariance[3,i] = cov(xval[is.q3],yval[is.q3])
        #q4
        is.q4 = (xval>(-hole_sizes[i]/xval) & yval<(-hole_sizes[i]/xval))
        occurrence[4,i] = sum(is.q4,na.rm=TRUE)
        product[4,i] = mean(xval[is.q4] * yval[is.q4],rm=TRUE)
        covariance[4,i] = cov(xval[is.q4],yval[is.q4])
    }
    #exuberance and organization ratio
    if (orient == "-") {
        exub=(product[1,]+product[3,])/(product[2,]+product[4,])
        or=(occurrence[1,]+occurrence[3,])/(occurrence[2,]+occurrence[4,])
    } else if (orient == "+") {
        exub=(product[2,]+product[4,])/(product[1,]+product[3,])
        or=(occurrence[2,]+occurrence[4,])/(occurrence[1,]+occurrence[3,])
    } else {
        warning("The orientation has to be either + or -.")
    }
    qa_out=list("hole_sizes"=hole_sizes,
            "occurrence"=occurrence,
            "product"=product,
            "covariance"=covariance,
            "covariance_total"=covariance_total,
            "correlation_total"=correlation_total,
            "exuberance"=exub,
            "organization_ratio"=or,
            "meta"="Output format: rows represent the quadrants Q1, Q2, Q3, Q4 -- columns represent selected hole sizes")
    if (plot==TRUE) {
        plot_quadrant_analysis(xval,yval,do_normalization=FALSE,...)
        print(qa_out)
    }
    return(qa_out)
}

#' Plotting Quadrant Analysis
#'
#'@description Plots two vectors in the framework of quadrant analysis with 2d kernel density estimation (optional)
#'@param xval values of x variable (vector)
#'@param yval values of y variable (vector)
#'@param do_normalization should the values be normalized? i.e. \code{(x-mean(x))/sd(x)}, default: \code{do_normalization=TRUE}
#'@param hole_sizes vector containing desired hole sizes (integers >= 0), default: \code{hole_sizes=c(1,2)}
#'@param plot_kde2d should the contour lines of the 2d kernel density estimation be plotted? default \code{plot_kde2d = TRUE}
#'@param contours vector containing levels of contour lines for 2d kernel density estimation, only used if \code{plot_kde2d = TRUE}, default: \code{contours=10^(-3:3)}
#'@param print_fit should the fit summary from the linear regression be printed? default: \code{print_fit=TRUE}
#'@param ... arguments passed to plot function
#'@return no return
#'@export
#'
#'@importFrom MASS kde2d
#' 
#'@examples
#'a=rnorm(100)
#'b=rnorm(100)
#'plot_quadrant_analysis(a,b)
#'
plot_quadrant_analysis=function(xval,yval,do_normalization=TRUE,hole_sizes=c(1,2,3),plot_kde2d=TRUE,contours=10^(-3:3),print_fit=TRUE,...) {
    if (do_normalization) {
        xval=(xval-mean(xval,na.rm=TRUE))/sd(xval)
        yval=(yval-mean(yval,na.rm=TRUE))/sd(yval)
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
    if (plot_kde2d==TRUE) {
        nc=length(contours)
        lab=colorRampPalette(c("blue3","red3"), space = "Lab")
        notna=(!is.na(xval) & !is.na(yval))
        kde=MASS::kde2d(xval[notna],yval[notna])
        contour(kde$x,kde$y,kde$z,levels=contours,col=lab(nc)[1:nc],add=TRUE,lwd=1,drawlabels=TRUE)
    }
    #draw holes
    xs=seq(-10,10,0.1)
    for (h in hole_sizes) {
        points(xs,h/xs,col=1,type="l",lty=2,lwd=1,pch=pch)
        points(xs,-h/xs,col=1,type="l",lty=2,lwd=1,pch=pch)
    }
    legend("bottomleft",legend=c("linear fit","hyperbolic holes"), col=c("darkred",1),lwd=c(2,1),lty=c(1,2))
}



#' Two-point Correlation Function: Tool to study spatial characteristics of coherent structures
#'
#'@description 
#'@param ts1 timeseries at point 1 (high-resolution)
#'@param ts2 timeseries at point 2 (high-resolution)
#'@return scalar giving the two-point correlation R(ts1,ts2)
#'@export
#' 
#'@examples
#'set.seed(5)
#'ts1=rnorm(100)
#'ts2=rnorm(100)
#'cor_2point=calc_2point_cor(ts1,ts2)
#'
calc_2point_cor=function(ts1,ts2) {
    scor=cor(ts1,ts2,use="pairwise.complete.obs")
    return(scor/(sd(ts1,na.rm=TRUE)*sd(ts2,na.rm=TRUE)))
}


