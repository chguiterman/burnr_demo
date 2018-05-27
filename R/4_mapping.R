## Demo script for R burnr package
## Chris Guiterman
## chguiterman@email.arizona.edu
## 2018-5-27

## Module 4: Mapping

library(burnr)

data("pgm")
data("pgm_meta")

# a couple useful libraries
library(rgdal)
library(maptools)

# Read in the shrubfield polygon for PGM
pgm_shrub <- readShapePoly('Data/PGM_shp/pgm_shrub_poly.shp', 
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
  #geom_polygon(data=pgm_shrub, aes(x=long, y=lat), fill=NA, colour='red') +
  scale_color_manual("Evidence type", values=col.man)

#' A nice example of broader-scale mapping with facet plotting is provided in
#' Malevich et al. 2018, Dendrochronologia
#' With code in the "bando_map" folder:
#' https://github.com/brews/burnr_2018_manuscript_figures



