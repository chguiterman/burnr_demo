## Demo script for R burnr package
## Chris Guiterman
## chguiterman@email.arizona.edu
## 2018-5-27

## Module 2: Graphics

library(burnr)
data(pgm)

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

#' This script will draw Figure 6 in,
#' Guiterman, C.H., E.Q. Margolis, C.D. Allen, D.A. Falk, 
#' and T.W. Swetnam. (2017). Long-term persistence and 
#' fire resilience of oak shrubfields in dry conifer 
#' forests of northern New Mexico. Ecosystems. 
#' https://doi.org/10.1007/s10021-017-0192-2

#' FHX2 files 
sns <- read_fhx('Data/SNS.fhx')
snn <- read_fhx('Data/SNN.fhx')
spc <- read_fhx('Data/SPC.fhx')
mpb <- read_fhx('Data/MPB.fhx')
rdc <- read_fhx('Data/RDC.fhx')

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

#' To add the composites from individual sites to that site within the facet, see the wiki page:
#' https://github.com/ltrr-arizona-edu/burnr/wiki/burnr-Cookbook#facetcomposite
