# Reddy (dev): Analyzing Turbulence Data
<!-- badges: start --> 
[![CRAN status](https://www.r-pkg.org/badges/version/Reddy)](https://cran.r-project.org/package=Reddy)
[![Last Commit](https://img.shields.io/github/last-commit/noctiluc3nt/Reddy)](https://github.com/noctiluc3nt/Reddy)
[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
<!-- badges: end -->

R-Package (under development) to automatize standard analyses of turbulence data (e.g., from eddy-covariance measurements). 


## How to use the package

### Installation
The package Reddy can be installed directly from github (current development version):
```
    devtools::install_git("https://github.com/noctiluc3nt/Reddy")
```

### Functionality
The Reddy package provides functions for the post-processing, analysis and evaluation of eddy-covariance measurements (see schematic), which are described in the [manual](https://github.com/noctiluc3nt/Reddy/tree/main/inst/figures/Reddy-manual.pdf) and are divided into the following scripts:
- `anisotropy.R`: invariant analysis of the Reynolds stress tensor, calculation of turbulence anisotropy and visualization in a barycentric map
- `auxilliary.R`: collection of some useful generic functions for the evaluation (e.g., discrete binning, cross-correlation maximization)
- `constants.R`: constants used for calculations (internal)
- `diagnostics-meteorology.R`: calculation of "background-meteorology" quantities (e.g., clear-sky index)
- `diagnostics-turbulence.R`: calculation of some standard turbulence diagnostics (e.g., friction velocity, TKE, turbulence intensity, stability parameter)
- `ec-processing.R`: collection of functions for post-processing and quality control of eddy-covariance measurements
- `flux-footprint.R`: calculation and visualization of 2D flux footprint (FFP, Kljun et al., 2015)
- `multiresolution-decomposition.R`: calculation and visualization of multiresolution decomposition (MRD, Vickers and Mahrt, 2003)
- `quadrant-analysis.R`: calculation and visualization of quadrant analysis
- `surface-energy-balance.R`: visualization of surface energy balance, residual flux and closure ratio

<image src="./inst/figures/schema.png">


### Usage
A detailed tutorial containing several jupyter notebooks showcasing the most important features of Reddy along with some theoretical background information can be found [here](https://github.com/noctiluc3nt/ec_analyze).<br><br>
**Plot examples:**<br>
<image src="./inst/figures/reddy_showcase.png">

### Workflow for package building
- create documentation (Rd files) with roxygen2: `roxygen2::roxygenize("./Reddy")` 
- build package: `devtools::build("./Reddy")` or via terminal `R CMD build Reddy`
- check package: `devtools::check("./Reddy")` or via terminal `R CMD check Reddy`

### Other packages for processing of eddy-covariance data
- [EddyPro®](https://github.com/LI-COR-Environmental/eddypro-engine): Post-processing of eddy-covariance data (LI-COR Biosciences)
- [REddyProc](https://cran.r-project.org/web/packages/REddyProc/index.html): Post-processing of (half-)hourly eddy-covariance measurements
- [openeddy](https://github.com/lsigut/openeddy): Post-processing of eddy-covariance data, aligned with REddyProc
- [RFlux](https://github.com/domvit81/RFlux/): GUI for post-processing of eddy-covariance raw data by calling [EddyPro®](https://github.com/LI-COR-Environmental/eddypro-engine) (LI-COR Biosciences, 2019)
- [ONEFlux](https://github.com/fluxnet/ONEFlux) ("Open Network-Enabled Flux processing pipeline"): Post-processing of (half-)hourly eddy-covariance data used to create the [FLUXNET2015 dataset](https://fluxnet.org/data/fluxnet2015-dataset/)
- [icoscp](https://pypi.org/project/icoscp/): accessing data from [ICOS](https://www.icos-cp.eu/) (Integrated Carbon Observing System) [data portal](https://data.icos-cp.eu/portal/#%7B%22filterCategories%22:%7B%22project%22:%5B%22icos%22%5D,%22level%22:%5B1,2%5D,%22stationclass%22:%5B%22ICOS%22%5D%7D%7D)
- something missing: please create an issue


## Literature
- Foken, T. (2017). Micrometeorology, Springer-Verlag Berlin Heidelberg, doi: 10.1007/978-3-642-25440-6.
- Kljun, N. and Calanca, P. and Rotach, M. W. and Schmid, H. P. (2015). A simple two-dimensional parameterisation for Flux Footprint
Prediction (FFP), Geoscie. Model Dev., 8, 3695-3713.
- Vickers, D. and Mahrt, L. (2003). The Cospectral Gap and Turbulent Flux Calculations. Journal of Atmospheric and Oceanic Technology, 20:660-672.
