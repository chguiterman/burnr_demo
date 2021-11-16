## Demo script for R burnr package
## Chris Guiterman
## chguiterman@email.arizona.edu
## 2018-5-27

## Module 4: Mapping

library(burnr)
library(sf)
library(dplyr)
library(ggplot2)

# Load datasets
data("pgm")
data("pgm_meta")

# Load a shapefile with the shrubfield boundary
pgm_shrub <- read_sf("Data/PGM_shp/pgm_shrub_poly.shp") 

# check it out via plot
p <- ggplot() +
  geom_sf(data = pgm_shrub, fill = "grey95")

p

# Combine the fhx object with the metadata so we can map the ring data
pgm_loc <- inner_join(pgm, pgm_meta, by = c("series" = "TreeID")) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326)

head(pgm_loc)

# Check out the map
p + geom_sf(data = pgm_loc)


## The objective of the PGM site collection was to test how trees recorded a
## high-severity fire event in 1993 (the Buchannon Fire). We collected 4 lines
## of tree-ring evidence: fire scars, pith dates on resprouting shrubs (Gambel
## oak), death dates on fire-killed trees, and growth changes on surviving trees
## that did not record a fire scar - typically these trees had lifted crowns.


pgm_buchannon <- pgm_loc %>% 
  filter(year %in% 1992:1993,
         ! rec_type %in% "recorder_year") %>% 
  mutate(evidence = case_when(
    rec_type == "pith_year" ~ "Pith date",
    rec_type == "bark_year" ~ "Death date",
    grepl("_fs", .$rec_type) ~ "Fire scar",
    grepl("_fi", .$rec_type) ~ "Growth change"
  ))

## We can map these evidences across the site for the years 1992 and 1993

p + geom_sf(data = pgm_buchannon,
            aes(color = evidence),
            size = 2) +
  scale_color_manual(values = c("red", "orange", "purple", "black"),
                     name = "Line of evidence") +
  ggtitle("Tree-ring evidence for 1993 type conversion") +
  theme_bw()

#' A nice example of broader-scale mapping with facet plotting is provided in
#' Malevich et al. 2018, Dendrochronologia
#' With code in the "bando_map" folder:
#' https://github.com/brews/burnr_2018_manuscript_figures



