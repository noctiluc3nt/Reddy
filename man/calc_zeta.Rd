% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics-turbulence.R
\name{calc_zeta}
\alias{calc_zeta}
\title{Stability Parameter}
\usage{
calc_zeta(z, L)
}
\arguments{
\item{z}{measurement height [m]}

\item{L}{Obukhov length [m] (e.g., from \code{calc_L})}
}
\value{
stability parameter [-]
}
\description{
Calculates dimensionless stability parameter from Obukhov length and measurement height, i.e. \code{zeta = z/L}
}
\examples{
calc_zeta(2,-1) #unstable
calc_zeta(2,1) #stable

}
