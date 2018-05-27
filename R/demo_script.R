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
series_stats(pgm)

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
snn <- read_fhx('Data/SNN.fhx')
spc <- read_fhx('Data/SPC.fhx')
mpb <- read_fhx('Data/MPB.fhx')
rdc <- read_fhx('Data/RDC.fhx')
#' Guiterman, C.H., E.Q. Margolis, C.D. Allen, D.A. Falk, 
#' and T.W. Swetnam. (2017). Long-term persistence and 
#' fire resilience of oak shrubfields in dry conifer 
#' forests of northern New Mexico. Ecosystems. 
#' https://doi.org/10.1007/s10021-017-0192-2


# Graphics ----------------------------------------------------------------

#' Make a simple plot
plot(pgm)

#' try a different sorting
pgm <- sort(pgm, decreasing = FALSE, sort_by = 'first_year')

plot(pgm, ylabels = FALSE)

#' add a composite rug to the plot
plot(pgm, ylabels=FALSE, composite_rug=TRUE)

#' add a legend
plot(pgm, ylabels=FALSE, composite_rug=TRUE, plot_legend = TRUE)

##' Start making it cool

#' Color by species using a separate metadata file
data("pgm_meta")
head(pgm_meta)

plot(pgm, color_group = pgm_meta$SpeciesID, color_id = pgm_meta$TreeID,
     plot_legend = TRUE, ylabels=FALSE)

#' far more functionality using ggplot options with plot_demograph() 
library(ggplot2)

p <- plot_demograph(pgm, color_group = pgm_meta$SpeciesID, color_id = pgm_meta$TreeID,
                    plot_legend = TRUE, ylabels=FALSE)
p

#' Change colors
levels(pgm_meta$SpeciesID)
cols <- c('purple', 'brown','red', 'blue', 'green4')

p + scale_color_manual(values = cols)

#' Annotate patterns

p + scale_color_manual(values = cols) +
  annotate('rect', xmin = 1871, xmax = 1993, ymin = 0, ymax = 45, alpha = .3,
           fill = 'blue2')

# Advanced graphics -------------------------------------------------------

#' This script will draw Figure 6 in Guiterman et al. (2017, Ecosystems)

#' Read in associated metadata
trees <- read.csv('Data/Shrubfield_meta.csv', row.names=1)

#' sort objects
mpb <- sort(mpb, decreasing=FALSE, sort_by = "first_year")
rdc <- sort(rdc, decreasing=FALSE, sort_by = "first_year")
snn <- sort(snn, decreasing=FALSE, sort_by = "first_year")
sns <- sort(sns, decreasing=FALSE, sort_by = "first_year")
spc <- sort(spc, decreasing=FALSE, sort_by = "first_year")

#' combine files
all.fhx <- mpb + rdc + snn + sns + spc

#' Identify sites for faceting and filter metadata
trees$SiteID <- substr(trees$series, start=1, stop=3)
trees <- trees[trees$series %in% unique(all.fhx$series), ]

#' Create site composites
mpb.comp <- composite(mpb, filter_min_rec = 2, filter_min_events = 2, filter_prop = .10, injury_event = TRUE, comp_name = "MPB")
rdc.comp <- composite(rdc, filter_min_rec = 2, filter_min_events = 2, filter_prop = .10, injury_event = TRUE, comp_name = "RDC")
snn.comp <- composite(snn, filter_min_rec = 2, filter_min_events = 2, filter_prop = .10, injury_event = TRUE, comp_name = "SNN")
sns.comp <- composite(sns, filter_min_rec = 2, filter_min_events = 2, filter_prop = .10, injury_event = TRUE, comp_name = "SNS")
spc.comp <- composite(spc, filter_min_rec = 2, filter_min_events = 2, filter_prop = .10, injury_event = TRUE, comp_name = "SPC")

#' combine composites
all.comp <- spc.comp + sns.comp + snn.comp + rdc.comp + mpb.comp

all <- all.fhx + all.comp

#' Add composite metadata to trees table
comps <- data.frame(series = c("MPB", "RDC", "SNN", "SNS", "SPC"),
                    SpeciesID = "COMP", SamplePurpose = NA, Latitude = NA,
                    Longitude = NA, Easting = NA, Northing = NA, Elevation = NA,
                    SiteID = "comp", Patch_loc = "comp")
trees.all <- rbind(trees, comps)


#' Organize metadata, now that "COMP" is a Species ID
trees.all$SpeciesID <- factor(trees.all$SpeciesID)
levels(trees.all$SpeciesID)
#' Assign colors for each species
col.sel <- c("#8B5A2B", "#4DAF4A", "#E41A1C", "#984EA3", "#405ECB", "black")
#' Organize the site levels for plotting
trees.all$SiteID <- factor(trees.all$SiteID, levels = c("MPB", "RDC", "SNN", "SNS", "SPC", "comp"))

#' Make the graph
p2 <- plot_demograph(all, color_group=trees.all$SpeciesID, color_id=trees.all$series,
                     facet_group=trees.all$SiteID, facet_id=trees.all$series,
                     event_size = c(1.75,1,1), #size of fs, inj, pith/bark
                     composite_rug=FALSE, rugdivide_pos = 0, #plot composite = YES and location
                     plot_legend = TRUE, ylabels=FALSE, #no treeID labels
                     yearlims= c(1500,1995),
                     facet_type='grid')
p2 + scale_color_manual(values=col.sel, breaks = c("JUSC", "PIPO", "PIST", "POTR", "QUGA")) + 
  guides(linetype = FALSE, shape = FALSE, size = FALSE) +
  theme(legend.justification = "center", legend.position = "top",
        legend.direction="horizontal", legend.background=element_rect(fill='white'),
        legend.margin=margin(0, unit="pt"),
        strip.text.y = element_text(size = 6)) +
  scale_x_continuous(sec.axis = dup_axis())
#' If you'd like to save the graph
ggsave("Output/Shrubfields_Fig6.tiff", device="tiff", width=4, height=6, dpi=300)


