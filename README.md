# Reddy (dev): Analyzing Eddy-Covariance Measurements

R-Package (under development) for analyzing and visualizing turbulent data, e.g., from eddy-covariance measurements.

## Scripts
- `anisotropy.R`: invariant analysis of the Reynolds stress tensor, calculation of turbulence anisotropy and visualization in a barycentric map
- `auxillary.R`: some (maybe) useful auxillary functions
- `mrd.R`: calculation and visualization of multiresolution decomposition 
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
- Vickers, D. and Mahrt, L. (2003). The Cospectral Gap and Turbulent Flux Calculations. Journal of Atmospheric and Oceanic Technology, 20:660-672.