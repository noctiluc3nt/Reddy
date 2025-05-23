% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-utils.R
\name{df_dy}
\alias{df_dy}
\title{df_dy}
\usage{
df_dy(fld, yres = 1)
}
\arguments{
\item{fld}{input field with dimension (x,y)}

\item{yres}{resolution in y-direction}
}
\value{
y-derivative of fld (same dimensions)
}
\description{
Calculates y-derivative for equidistant grid
}
\examples{
set.seed(5)
field=matrix(rnorm(16),ncol=4)
df_dy(field,10)

}
