#' Plottting of surface energy balance and calculate surface energy balance unclosure
#'
#'@description Plottting of surface energy balance and calculate surface energy balance unclosure as residual flux and closure ratio
#'@param sw_out outcomin shortwave radiation [W/m^2] (as vector of time)
#'@param lw_in incoming longwave radiation [W/m^2] (as vector of time)
#'@param lw_out outgoing longwave radiation [W/m^2] (as vector of time)
#'@param sh sensible heat flux [W/m^2] (as vector of time) -- if measured
#'@param lh latent heat flux [W/m^2] (as vector of time) -- if measured
#'@param gh ground heat flux [W/m^2] (as vector of time) -- if measured
#'@param time_vector times used as x-axis labels (optional)
#'@param print_fit should the fit summary be printed? default: TRUE
#'
#'@return no return
#'@export
#'
#'@examples
#'
plot_seb = function(sw_in,sw_out,lw_in,lw_out,sh=NULL,lh=NULL,gh=NULL,time_vector=NULL,print_fit=TRUE,...) {
    #check which fluxes are given
    if (is.null(sh)) sh = 0
    if (is.null(lh)) lh = 0
    if (is.null(gh)) gh = 0
    #balances
    rad_balance = sw_in-sw_out+lw_in-lw_out #radiation balance
    seb_balance = rad_balance-gh-sh-lh #surface energy balance (residual flux)
    cr = (rad_balance-gh)/(sh+lh) #closure ratio
    #plot
    n=length(rad_balance)
    if (is.null(time_vector)) time_vector = 1:n
    if (!exists("ylim")) ylim=c(-100,500)
    if (!exists("pch")) pch=20
    if (!exists("cex")) cex=2
    plot(time_vector,rad_balance,type="l",lwd=2,col="gray50",main="Surface Energy Balance",ylim=ylim,ylab="surface energy balance fluxes",xlab="index",...)
    points(time_vector,seb_balance,type="l",lwd=3,col=1)
    points(time_vector,sh,type="l",lwd=2,col="orangered")
    points(time_vector,lh,type="l",lwd=2,col="blue3")
    points(time_vector,gh,type="l",lwd=2,col="brown")
    legend("topright",legend=c("R","GH","SH","LH","Res"),col=c("gray50","brown","orangered","blue3",1),lwd=c(2,2,2,2,3),lty=1,bg="white")
    plot(rad_balance-gh,sh+lh,col=rgb(0,0,0,0.4),pch=pch,cex=cex,main="Surface Energy Balance Closure",xlab="R-G",ylab="SH+LH",...)
    y=sh+lh
    x=rad_balance-gh
    fit=lm(y ~ x)
    if (print_fit==TRUE) print(summary(fit))
    abline(fit,lwd=2,col=2,lty=2)
    abline(0,1,lwd=2,lty=1)
    grid()
}