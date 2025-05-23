% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-utils.R
\name{df_dx}
\alias{df_dx}
\title{df_dx}
\usage{
df_dx(fld, xres = 1)
}
\arguments{
\item{fld}{input field with dimension (x,y)}

\item{xres}{resolution in x-direction}
}
\value{
x-derivative of fld (same dimensions)
}
\description{
Calculates x-derivative for equidistant grid
}
\examples{
set.seed(5)
field=matrix(rnorm(16),ncol=4)
df_dx(field,10)

}
