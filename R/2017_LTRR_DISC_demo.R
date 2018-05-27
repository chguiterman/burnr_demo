## A demo of the burnr library in R
## C. Guiterman 05-31-17
## chguiterman@email.arizona.edu


# Install burnr -----------------------------------------------------------

# Download the library from cran
install.packages("burnr")

# load the library
library(burnr)


# Load and write fire history files --------------------------------------

# fire history data is often stored in FHX text files. These can read directly into R using 
# the read_fhx function
example <- read_fhx('example.fhx')

# Can be done directly from a web server, such as the IMPD
tt <- read_fhx('https://www1.ncdc.noaa.gov/pub/data/paleo/firehistory/firescar/northamerica/usber001.fhx')

# load some data from within the package - we have two datsets available
data("pgm")
data("lgr2")

# The data are stored in a simple data frame
head(pgm)

# There are 19 record types for each year
levels(pgm$rec_type)

series_names(pgm)

# you can combine the two datasets
all <- pgm + lgr2

# write out your new fhx file in FHX format
write_fhx(all, "newFile.fhx")


# Graphics ----------------------------------------------------------------

# make a simple plot
plot(pgm)

# try a different sorting
pgm <- sort(pgm, decreasing = FALSE, sort_by = 'first_year')

plot(pgm, ylabels = FALSE)

# add a composite rug to the plot
plot(pgm, ylabels=FALSE, composite_rug=TRUE)

# add a legend
plot(pgm, ylabels=FALSE, composite_rug=TRUE, plot_legend = TRUE)

## Start making it cool

# Color by species using a separate metadata file
data("pgm_meta")
head(pgm_meta)

plot(pgm, color_group = pgm_meta$SpeciesID, color_id = pgm_meta$TreeID,
     plot_legend = TRUE, ylabels=FALSE)

# far more functionality using ggplot options with plot_demograph 
library(ggplot2)

p <- plot_demograph(pgm, color_group = pgm_meta$SpeciesID, color_id = pgm_meta$TreeID,
                    plot_legend = TRUE, ylabels=FALSE)
p

# Change colors
levels(pgm_meta$SpeciesID)
cols <- c('purple', 'brown','red', 'blue', 'green4')

p + scale_color_manual(values = cols)

# Annotate patterns

p + scale_color_manual(values = cols) +
  annotate('rect', xmin = 1871, xmax = 1993, ymin = 0, ymax = 45, alpha = .3,
           fill = 'blue2')

# plot multiple sites or plots together using facets

## Read in data from Jemez shrubfield sites. Data stored by R from save.image()

load('C:/Users/chris.guiterman/Dropbox/Teaching/2017_Presession/Shrubfields.RData')

p2 <- plot_demograph(all, color_group=trees.all$SpeciesID, color_id=trees.all$series,
                     facet_group=trees.all$SiteID, facet_id=trees.all$series,
                     event_size = c(1.75,1,1), #size of fs, inj, pith/bark
                     composite_rug=FALSE, rugdivide_pos = 0, #plot composite = YES and location
                     plot_legend = TRUE, ylabels=FALSE, #no treeID labels
                     yearlims= c(1500,1995),
                     facet_type='grid')
p2 + scale_color_manual(values=col.sel, breaks = c("JUSC", "PIPO", "PIST", "POTR", "QUGA")) + 
  guides(linetype = FALSE, shape = FALSE, size = FALSE) +
  theme(legend.justification = "center", legend.position = c(0.073, 0.5), 
        legend.direction="vertical", legend.background=element_rect(fill='white'), 
        legend.box="vertical")


# Analyses ----------------------------------------------------------------

# See some basic stats for each tree (series)
head(series_stats(pgm))

# Create a composite
pgm.comp <- composite(pgm, filter_min_events = 1, filter_prop = 0, filter_min_rec = 1,
                      comp_name = 'PGM')
plot(pgm.comp)

## Intervals
intervals(pgm.comp)

intervals(pgm.comp, densfun = "lognormal")

# Make it an intervals object and plot the data
plot(intervals(pgm.comp), binwidth = 5)

# Check out seasonality
count_event_position(pgm)

## Superposed Epoch Analysis

# Read in PDSI timeseries for pgm
data("pgm_pdsi")

head(pgm_pdsi)

(pgm.sea <- sea(x = pgm_pdsi, event = pgm.comp, nbefore = 4, nafter = 2))

str(pgm.sea)

plot(pgm.sea)

# Replicate the SEA from Dr. Falk's demo:
zmt <- read_fhx('http://www.ltrr.arizona.edu/~sheppard/presession/ZMT.FHX')

pdsi <- read.csv('http://www.ltrr.arizona.edu/~sheppard/presession/NADA GP119 1350-1900.csv', 
                 row.names = 1)

zmt.comp <- composite(zmt, filter_prop = .25, filter_min_rec = 1, filter_min_events = 2)

plot(sea(pdsi, zmt.comp))

# Extensions in the R ecosystem -------------------------------------------

## Mapping

library(rgdal)
library(maptools)

# Read in the shrubfield polygon for PGM
pgm_shrub <- readShapePoly('~/burnr_develop/PGM_shp/pgm_shrub_poly.shp', 
                           proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))

plot(pgm_shrub)

# Merge PGM metadata and fhx objects
pgm.loc <- merge(pgm, pgm_meta, by.x = "series", by.y = "TreeID")

head(pgm.loc)

# Designate pgm.loc at spatial data object
coordinates(pgm.loc) <- c("Longitude", "Latitude")
# Specify the coordinate projection
proj4string(pgm.loc) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84')

# plot tree locations
plot(pgm_shrub)
points(pgm.loc)

# Add column to specify colors for different rec_types
pgm.loc$col <- 'grey80' # will leave recording status as grey later
pgm.loc[grep('_fs', pgm.loc$rec_type), 'col'] <- 'red' # Make fire scars red
pgm.loc[grep('_fi', pgm.loc$rec_type), 'col'] <- 'orange' # Make injuries orange
pgm.loc[grep('pith_year', pgm.loc$rec_type), 'col'] <- 'green' # Make pith dates green
pgm.loc[grep('bark_year', pgm.loc$rec_type), 'col'] <- 'black' # Make bark dates black

## Draw map for the effects high-severity burning in spring 1993
plot(pgm_shrub)
# add death dates dating to 1992 (last complete ring)
points(pgm.loc[pgm.loc$year == 1992, ], col=pgm.loc[pgm.loc$year == 1992, ]$col, pch=20)
# add fire scars, growth chnges, and  pith dates
points(pgm.loc[pgm.loc$year == 1993, ], col=pgm.loc[pgm.loc$year == 1993, ]$col, pch=20)
title('Tree-ring evidence for 1993 type conversion')
legend('bottomleft', c('Pith date', 'Fire scar', 'Growth change', 'Death date'),
       col=c('green', 'red', 'orange', 'black'), pch=20)

## Use a different background
library(ggmap)
library(ggplot2)


map <- get_map(location = c(lon = mean(pgm_meta$Longitude), lat = mean(pgm_meta$Latitude)), 
               zoom = 16, source = "google", maptype = "satellite", crop = FALSE)
p <- ggmap(map)
p + geom_point(aes(x = Longitude, y = Latitude), data = pgm_meta, alpha = 1, 
               color = "darkred", size = 5) + theme_bw()

# Highlight the 1993 Buchannon fire
pgm2 <- data.frame(pgm.loc)
pgm2[grep('_fs', pgm2$rec_type), 'EvidenceType'] <- 'Fire scar'
pgm2[grep('_fi', pgm2$rec_type), 'EvidenceType'] <- 'Growth change'
pgm2[grep('pith_year', pgm2$rec_type), 'EvidenceType'] <- 'Pith date'
pgm2[grep('bark_year', pgm2$rec_type), 'EvidenceType'] <- 'Death date'
pgm2$EvidenceType <- factor(pgm2$EvidenceType)
col.man <- c('black', 'red', 'orange', 'green')
# Subset the data to include only the evidences that align with the 1993 Buchanan Fire
pgm2.buchanan <- pgm2[pgm2$year == 1992 | pgm2$year == 1993, ]


p + geom_jitter(data=pgm2.buchanan, aes(x=Longitude, y=Latitude, col=EvidenceType), 
                shape=19, size=3, width=.0005, height=.0005) +
  geom_polygon(data=pgm_shrub, aes(x=long, y=lat), fill=NA, colour='red') +
  scale_color_manual("Evidence type", values=col.man)