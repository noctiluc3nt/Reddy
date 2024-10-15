#!/usr/bin/Rscript

roxygen2::roxygenize("../Reddy")

devtools::build("../Reddy")
devtools::check("../Reddy")
