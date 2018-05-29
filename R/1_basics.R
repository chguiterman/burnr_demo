## Demo script for R burnr package
## Chris Guiterman
## chguiterman@email.arizona.edu
## 2018-5-27

## Module 1: Basics

# Install burnr -----------------------------------------------------------

# Download the library from cran
install.packages("burnr")

# load the library
library(burnr)

#' If you're looking for new features, check out the development version of burnr before it goes up to CRAN

library(devtools)
install_github("ltrr-arizona-edu/burnr@dev")


# Obtaining and loading data ----------------------------------------------

#' burnr is pre-installed with two fire histories, each from the Jemez Mountains in New Mexico. 
#' Once in the R environment, these datasets are known as fhx objects. To runt he burnr functions,
#' the fire history datasets must be fhx objects.

data("lgr2")
#' Margolis, E. Q., Malevich, S. B., 2016. Historical dominance of low-severity
#' fire in dry and wet mixed-conifer forest habitats of the endangered terrestrial 
#' Jemez Mountains salamander (Plethodon neomexicanus). Forest Ecology and Management 375, 12-26.

data("pgm")
#' Guiterman, C. H., Margolis, E. Q., Swetnam, T. W., 2015. Dendroecological
#' Methods For Reconstructing High-Severity Fire In Pine-Oak Forests. Tree-
#' Ring Research 71 (2), 67-77.

#' These data can be viewed in several ways.
View(pgm)
head(pgm)

#' tree-level summaries are provided by
head(series_stats(pgm))

#' series (tree) names are provided by
series_names(pgm)

#' FHX objects are rather simple, with 3 columns for "year", "series", and "rec_type"
#' The "rec_type" is a catch-all for the tree-ring attribute for a given year/series
#' There are currently 19 rec_type levels
levels(pgm$rec_type)

#' Oportuities to modify and edit fhx objects abound.
delete(lgr2, s = 'LGR46')             # Remove series.
delete(lgr2, yr = 1752)               # Remove year from all series.
delete(lgr2, s = 'LGR46', yr = 1752)  # Remove year from select series.
#' Use help() to find more information on these functions. 
#' Be sure to check out get_series(), get_year(), and sort.fhx().

#' Note that none of these actions are saving to the lgr2 object. To save
#' changes, make a new object, or re-assign the existing object

#' burnr is also capable of reading and writing FHX2 file formats.
#' These can obtained from the International Multiproxy Paleofire Databank.
url <- "https://www1.ncdc.noaa.gov/pub/data/paleo/firehistory/firescar/northamerica/"
pmr <- read_fhx(paste0(url, 'uspmr001.fhx'))
pme <- read_fhx(paste0(url, 'uspme001.fhx'))
pmw <- read_fhx(paste0(url, 'uspmw001.fhx'))

#' Combine all three sites into a single fhx object, and save the new FHX2 file
pm_all <- pmr + pme + pmw
write_fhx(pm_all, "Output/PajaritoMountain_allSites.fhx")

#' FHX2 files stored on the hard drive:
sns <- read_fhx('Data/SNS.fhx')

#' Guiterman, C.H., E.Q. Margolis, C.D. Allen, D.A. Falk, 
#' and T.W. Swetnam. (2017). Long-term persistence and 
#' fire resilience of oak shrubfields in dry conifer 
#' forests of northern New Mexico. Ecosystems. 
#' https://doi.org/10.1007/s10021-017-0192-2
