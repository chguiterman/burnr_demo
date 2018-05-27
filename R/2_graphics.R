## Demo script for R burnr package
## Chris Guiterman
## chguiterman@email.arizona.edu
## 2018-5-27

## Module 2: Graphics

library(burnr)


# Basic plotting ----------------------------------------------------------

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
trees <- trees[trees$series %in% series_names(all.fhx), ]

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


# FHAES-style graph -------------------------------------------------------

library(patchwork)
#' patchwork is available from
#' https://github.com/thomasp85/patchwork
#' devtools::install_github("thomasp85/patchwork")
#' It's still under development, and loading the library can
#' be a real pain. Be patient, and use iteration. It takes some
#' time, but there are some nice benefits

library(ggrepel) # For labeling
library(scales)

#' For data, let's use a single site from the shrubfields
snn_fhx <- read_fhx('Data/SNN.fhx')
snn_fhx <- sort(snn_fhx)

#' Start with percent trees scarred timeseries

snn_perc <- percent_scarred(snn_fhx)

#' We're going to add year labels to identify certain fire events, based on filtering
yr_labs <- snn_perc[snn_perc$num_scars > 1 & snn_perc$percent_scarred >= 25, ] 
wide_labs <- snn_perc[snn_perc$num_scars > 1 & snn_perc$percent_scarred >= 45, ] 

#' Patchwork operates by adding saved ggplot graph objects together, so we make each of the 3 sections
#' separatey

p <- ggplot() + 
  geom_col(data=snn_perc, aes(x=year, y=percent_scarred)) +
  xlim(1550, 2020) +
  geom_text(data=wide_labs, aes(x=year, y=percent_scarred+5, label=year, angle=35),
            size=2, nudge_x = 5) +
  scale_y_continuous(name = "% trees\nscarred", limits=c(0, 100), expand=c(0, 0)) +
  scale_x_continuous(position = "top", limits=c(1550, 2020)) +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
f <- plot_demograph(snn_fhx, yearlims = c(1550, 2020), composite_rug = TRUE,
                    filter_prop = 0.25, filter_min_events = 2, ylabels = FALSE,
                    plot_legend = TRUE) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = c(.15, .75))

l <- ggplot() +
  geom_text_repel(data=yr_labs, aes(x=year, 
                                    y=rep(1, nrow(yr_labs)),
                                    label = year, angle = -85), size=2.5,
                  direction = "y", segment.size = 0.5, segment.alpha=.5,
                  force=10) +
  scale_x_continuous(name="Year", limits = c(1550, 2020), breaks = seq(1500, 2000, 100)) +
  ylab("") + xlab("Year") +
  scale_y_continuous(limits=c(0, 1), expand=c(0, 0), breaks = c(0, 1), labels = NULL, 
                     minor_breaks = NULL, name = "\n\nComposite\nfire years") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.background = element_blank(),
        axis.ticks.y = element_blank(), panel.grid.major.x = element_line(linetype=1, color="grey90", size=.5),
        axis.line.x.bottom = element_line())

#' Stack them up, with adjustable heights

p + f + l + plot_layout(ncol=1, heights=c(1, 6, 1))

#' To save this, ggsave won't work, so go old-school

tiff('Output/FHAES-style.tiff', width=6, height=8, units='in', res=150, type='cairo')
p + f + l + plot_layout(ncol=1, heights=c(1, 6, 1))
dev.off()