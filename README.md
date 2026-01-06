# Reddy: A toolbox for analyzing eddy-covariance measurements
<!-- badges: start 
[![CRAN status](https://www.r-pkg.org/badges/version/Reddy)](https://cran.r-project.org/package=Reddy) -->
[![doc](https://img.shields.io/badge/Reddy-documentation-blue)](https://noctiluc3nt.github.io/ec_analyze/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Last Commit](https://img.shields.io/github/last-commit/noctiluc3nt/Reddy)](https://github.com/noctiluc3nt/Reddy)
<!--[![maintained](https://img.shields.io/badge/maintained-black)
<!-- badges: end -->


## Installation
The package Reddy can be installed directly from github:
```
    devtools::install_git("https://github.com/noctiluc3nt/Reddy")
```

## Documentation, Examples and Usage
The [manual](https://github.com/noctiluc3nt/Reddy/tree/main/inst/manual/Reddy-manual.pdf) describes all functions and the [gitbook](https://noctiluc3nt.github.io/ec_analyze/) provides background information and examples how to execute the core functions of this package. The examples are also provided as [jupyter notebooks](https://github.com/noctiluc3nt/ec_analyze/tree/main/notebooks) (based on [example data](https://github.com/noctiluc3nt/ec_analyze/tree/main/data)).

## Functionality and Scripts
The Reddy package provides functions for the post-processing, analysis and evaluation of eddy-covariance measurements, which are described in the [manual](https://github.com/noctiluc3nt/Reddy/tree/main/inst/manual/Reddy-manual.pdf) and are divided into the following scripts:
- `anisotropy.R`: invariant analysis of the Reynolds stress tensor, calculation of turbulence anisotropy and visualization in a barycentric map
- `auxilliary.R`: collection of some useful generic functions for the evaluation (e.g., discrete binning, cross-correlation maximization)
- `constants.R`: constants used for calculations (internal)
- `diagnostics-meteorology.R`: calculation of "background-meteorology" quantities (e.g., clear-sky index, vapor pressure deficit, unit conversions of moisture variables)
- `diagnostics-topography.R`: calculation of topographic and surface properties (e.g. topographic position index, slope-based parameter, roughness length)
- `diagnostics-turbulence.R`: calculation of some standard turbulence diagnostics (e.g., friction velocity, TKE, turbulence intensity, stability parameter)
- `ec-processing.R`: collection of functions for post-processing and quality control of eddy-covariance measurements
- `ec-processing-routine.R`: an example for a post-processing and quality control routine of eddy-covariance measurements utilizing the functions in `ec-processing.R` with storing the final averaged output data in an external file
- `flux-footprint.R`: calculation and visualization of 1D and 2D flux footprint based on the parameterization from Kormann & Meixner, 2001 and Kljun et al., 2015
- `model-utils.R`: collection of function to post-process NWP model output (e.g. flux deaccumulation, conversion of sigma levels to physical height)
- `multiresolution-decomposition.R`: calculation and visualization of multiresolution decomposition (MRD) and derivation of a suggestion for a suitable averaging time following the MRD-based algorithm from Vickers & Mahrt, 2003
- `ogive.R`: calculation and visualization of ogives, i.e. cumulative distribution functions
- `parameterizations-of-turbulence.R`: collection of functions used to calculated bulk closures, flux-profile and flux-variance relations (e.g. Richardson number, eddy viscosity, eddy conductivity, scaling functions) used in turbulence parameterizations in numerical weather prediction and climate models
- `quadrant-analysis.R`: calculation and visualization of quadrant analysis to study coherent structures and their organization
- `spectrum.R`: calculation and visualization of frequency (temporal) and wavenumber (spatial) spectra, possibility to bin them to compare the dissipation rates with theoretically expected slopes
- `structure-functions.R`: calculation of structure functions and autocorrelation
- `surface-energy-balance.R`: visualization of surface energy balance, residual flux and closure ratio

<!--
<image src="./inst/figures/schema2.png">



### Workflow for package building
- create documentation (Rd files) with roxygen2: `roxygen2::roxygenize("./Reddy")` 
- build package: `devtools::build("./Reddy")` or via terminal `R CMD build Reddy`
- check package: `devtools::check("./Reddy")` or via terminal `R CMD check Reddy`
-->

**Any issues or comments?** Create an issue [here](https://github.com/noctiluc3nt/Reddy/issues).

Contact: laura.mack@geo.uio.no

## Literature
- Mack, L., Berntsen, T.K., Vercauteren, N., Pirk, N. (2024). Transfer Efficiency and Organization in Turbulent Transport over Alpine Tundra. Boundary-Layer Meteorology 190, 38. DOI: https://doi.org/10.1007/s10546-024-00879-5



