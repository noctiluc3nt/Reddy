# Reddy (dev): Analyzing Turbulence Data
<!-- badges: start --> 
[![CRAN status](https://www.r-pkg.org/badges/version/Reddy)](https://cran.r-project.org/package=Reddy)
[![Last Commit](https://img.shields.io/github/last-commit/noctiluc3nt/Reddy)](https://github.com/noctiluc3nt/Reddy)
[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
<!-- badges: end -->

R-Package (under development) for analyzing turbulence data, e.g., from eddy-covariance measurements. 

## Scripts
contains currently
- `anisotropy.R`: invariant analysis of the Reynolds stress tensor, calculation of turbulence anisotropy and visualization in a barycentric map
- `constants.R`: constant used for calculations
- `ec-processing.R`: function for postprocessing of eddy-covariance measurements
- `footprint.R`: calculation and visualization of 2D flux footprint (FFP, Kljun et al., 2015)
- `multiresolution-decomposition.R`: calculation and visualization of multiresolution decomposition (MRD, Vickers and Mahrt, 2003)
- `quadrant-analysis.R`: calculation and visualization of quadrant analysis
- `turbulence-quantities.R`: calculation of some standard turbulence quantities

## Workflow
- create documentation (Rd files) with roxygen2: `roxygen2::roxygenize("./Reddy")` 
- build package: `devtools::build("./Reddy")` or via terminal `R CMD build Reddy`
- check package: `devtools::check("./Reddy")` or via terminal `R CMD check Reddy`

## Literature
- Kljun, N. and Calanca, P. and Rotach, M. W. and Schmid, H. P. (2015). A simple two-dimensional parameterisation for Flux Footprint
Prediction (FFP), Geoscie. Model Dev., 8, 3695-3713.
- Vickers, D. and Mahrt, L. (2003). The Cospectral Gap and Turbulent Flux Calculations. Journal of Atmospheric and Oceanic Technology, 20:660-672.

## other useful R-Packages for analyzing eddy-covariance data
- REddyProc: Post Processing of (Half-)Hourly Eddy-Covariance Measurements
