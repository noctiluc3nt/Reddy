#' Structure functions
#'
#'@description Calculates the structure function of given timeseries and given order, S := < (ts[t+dt], ts[t])^order >
#'@param ts timeseries
#'@param order order of the structure function, typically d = 2, 3, 4
#'@return structure function
#'@export
#'
#'@examples
#'ts=rnorm(100)
#'S2_ts=calc_structure_function(ts,2)
#'
calc_structure_function = function(ts,order=2) {
    return(mean(diff(ts)^order,na.rm=TRUE))
}


#' Two-point Correlation Function: Tool to study spatial characteristics of coherent structures
#'
#'@description Calculates two-poin correlation
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