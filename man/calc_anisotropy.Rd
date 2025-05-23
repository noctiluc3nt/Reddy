% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anisotropy.R
\name{calc_anisotropy}
\alias{calc_anisotropy}
\title{Invariant analysis of Reynolds stress tensor}
\usage{
calc_anisotropy(a11, a12, a13, a22, a23, a33, plot = FALSE)
}
\arguments{
\item{a11}{R11 element of Reynolds stress tensor: \code{u_sd^2}   (scalar or vector)}

\item{a12}{R12 element of Reynolds stress tensor: \code{cov(u,v)} (scalar or vector)}

\item{a13}{R13 element of Reynolds stress tensor: \code{cov(u,w)} (scalar or vector)}

\item{a22}{R22 element of Reynolds stress tensor: \code{v_sd^2}   (scalar or vector)}

\item{a23}{R23 element of Reynolds stress tensor: \code{cov(v,w)} (scalar or vector)}

\item{a33}{R33 element of Reynolds stress tensor: \code{w_sd^2}   (scalar or vector)}

\item{plot}{should the barycentric map be plotted? default \code{plot=FALSE}}
}
\value{
list containing \code{xb}, \code{yb}, \code{eta}, \code{xi}, all eigenvalues and eigenvectors (\code{eta}, \code{xi} are the coordinates of the Lumley triangle and \code{xb}, \code{yb} the coordinates of the barycentric map)
}
\description{
Invariant analysis of Reynolds stress tensor, calculation of Lumley and barycentric map coordinates and anisotropy
}
\examples{
calc_anisotropy(1,0,0,1,0,1) #isotropic
calc_anisotropy(1,0,1,1,0,1) #some anisotropy

}
