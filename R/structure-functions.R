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