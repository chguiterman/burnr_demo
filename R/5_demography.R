## Demo script for R burnr package
## Chris Guiterman
## chguiterman@email.arizona.edu
## 2018-5-27

## Module 5: Demography

library(burnr)
library(dplyr)
library(ggplot2)

#' We're going to use the shrubfield age-structure data

sns <- read_fhx('Data/SNSall.fhx')
snn <- read_fhx('Data/SNNall.fhx')
spc <- read_fhx('Data/SPCall.fhx')
mpb <- read_fhx('Data/MPBall.fhx')
rdc <- read_fhx('Data/RDCall.fhx')

all_fhx <- mpb + rdc + snn + sns + spc

#' Read in associated metadata
trees <- read.csv('Data/Shrubfield_meta.csv', row.names=1)

#' Identify sites in metadata
trees$SiteID <- substr(trees$series, start=1, stop=3)

#' Filter for age-structure trees
trees <- trees[trees$SamplePurpose == 'AGE', ]
age_fhx <- get_series(all_fhx, as.character(trees$series))

#' Filter for pith dates
age_fhx <- age_fhx[age_fhx$rec_type == "pith_year", ]

#' Merge tree-ring data and metadata
age_dat <- merge(age_fhx, trees, by='series')

#' Compile pith dates by site and species

site_dem <- age_dat %>% group_by(SiteID, SpeciesID, year) %>% 
  summarize(freq = n_distinct(series))

ggplot(site_dem) + geom_col(aes(x=as.numeric(year), y=freq, fill=SpeciesID)) +
  xlim(1830, 1950) +
  labs(x="Year", y="Number of stems") +
  facet_grid(SiteID~., scales='free_y')


#' Option to apply binning to years

bins <- seq(1450, 2020, 5)
age_dat$bin <- cut(age_dat$year, bins, right=FALSE)
age_dat$bin_yr <- bins[age_dat$bin]

site_dem_bin <- age_dat %>% group_by(SiteID, SpeciesID, bin_yr) %>% 
  summarize(freq = n_distinct(series))

ggplot(site_dem_bin) + geom_col(aes(x=bin_yr, y=freq, fill=SpeciesID)) +
  xlim(1830, 1950) +
  labs(x="Year", y="Number of stems") +
  facet_grid(SiteID~., scales='free_y')


