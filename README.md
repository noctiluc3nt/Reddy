# Reddy (dev): Analyzing Turbulence Data
<!-- badges: start --> 
[![CRAN status](https://www.r-pkg.org/badges/version/meteoEVT)](https://cran.r-project.org/package=meteoEVT)
[![Downloads (total)](http://cranlogs.r-pkg.org/badges/grand-total/meteoEVT?color=brightgreen)](https://cran.r-project.org/package=meteoEVT)
[![Last Commit](https://img.shields.io/github/last-commit/noctiluc3nt/meteoEVT)](https://github.com/noctiluc3nt/meteoEVT)
[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
<!-- badges: end -->

R-Package (under development) for analyzing turbulence data, e.g., from eddy-covariance measurements. 

## Scripts
- `anisotropy.R`: invariant analysis of the Reynolds stress tensor, calculation of turbulence anisotropy and visualization in a barycentric map
- `auxillary.R`: some (maybe) useful auxillary functions
- `footprint.R`: calculation and visualization of 2D flux footprint (FFP, Kljun et al., 2015)
- `mrd.R`: calculation and visualization of multiresolution decomposition (MRD, Vickers and Mahrt, 2003)
- `quantities_turbulence.R`
- `quantities_meteorology.R`
- `quadrant_analysis.R`
- `surface_energy_balance.R`
- scaling function?
- postprocessing?

## Workflow
- create documentation (Rd files) with roxygen2: `roxygen2::roxygenize("./Reddy")` 
- build package: `devtools::build("./Reddy")` or via terminal `R CMD build Reddy`
- check package: `devtools::check("./Reddy")` or via terminal `R CMD check Reddy`

## Literature
- Kljun, N. and Calanca, P. and Rotach, M. W. and Schmid, H. P. (2015). A simple two-dimensional parameterisation for Flux Footprint
Prediction (FFP), Geoscie. Model Dev., 8, 3695-3713.
- Vickers, D. and Mahrt, L. (2003). The Cospectral Gap and Turbulent Flux Calculations. Journal of Atmospheric and Oceanic Technology, 20:660-672.
