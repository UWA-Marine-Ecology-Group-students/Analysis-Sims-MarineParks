###
# Project: Hayley ecomix paper
# Data:    BRUVS
# Task:    Make publication ready plots
# author:  Kingsley Griffin, Hayley Sims
# date:    June 2022
##

library(sf)
library(raster)
library(reshape2)
library(ggplot2)
library(viridis)
library(patchwork)
library(ggnewscale)

# get predicted archetype rasters
predictions  <- list.files("output", "*predicted", full.names = TRUE)
pred_rasters <- stack(predictions)
plot(pred_rasters)

# get aus coastline
aus <- st_read("data/shapefiles/cstauscd_r.mif")
aus <- aus[aus$FEAT_CODE == "mainland", ]

# get NPZ and state sanctuary boundaries and subset to just our npz
amp_npz   <- st_read("data/shapefiles/AustraliaNetworkMarineParks.shp")
amp_npz   <- amp_npz[row.names(amp_npz) == 138, ]
wamp      <- st_read("data/shapefiles/test1.shp")
wamp      <- wamp[wamp$Name == "Cape Freycinet Sanctuary Zone", ]
wamp$Zone <- "Sanctuary Zone (WA)"

# set spatial reference system for aus outline and abrolhos to wgs84
st_crs(aus)  <- st_crs(amp_npz)
st_crs(wamp) <- st_crs(amp_npz)

# prepare predictions in long data frame for ggplot
pred_df            <- as.data.frame(pred_rasters, xy = TRUE, na.rm = TRUE)
colnames(pred_df)  <- gsub("predicted_", "", colnames(pred_df))                 # removing 'predicted'
pred_long          <- melt(pred_df, id.vars = c(1, 2))
head(pred_df)
head(pred_long)

# develop plot panes for each archetype
# set up various theme related things for plotting
thetheme <- theme_minimal() + 
  theme(axis.title = element_blank())                                           # theme elements
plotlims <- coord_sf(xlim = c(114.7, 115), 
                     ylim = c(-34.16, -34.0))                                   # crop area
amp_cols <- scale_colour_manual(values = c("National Park Zone" = "#7bbc63"))   # zone colours
wamp_col <- scale_colour_manual(values = c("Sanctuary Zone (WA)" = "palegreen"))     # state colours

# this is slow because rasters are large and detailed
p1 <- ggplot() +
  geom_raster(data = pred_df, aes(x, y, fill = Archetype1)) +
  scale_fill_viridis() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.2) +
  geom_sf(data = amp_npz, aes(colour = ZoneName), alpha = 4/5, fill = NA) +
  amp_cols +
  labs(colour = NULL) +
  new_scale_colour() +
  geom_sf(data = wamp, aes(colour = Zone), alpha = 4/5, fill = NA) +
  wamp_col +
  labs(colour = NULL) +
  guides(fill = guide_legend(order = 1)) +
  plotlims +
  thetheme
p1

p2 <- ggplot() +
  geom_raster(data = pred_df, aes(x, y, fill = Archetype2)) +
  scale_fill_viridis() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.2) +
  geom_sf(data = amp_npz, aes(colour = ZoneName), alpha = 4/5, fill = NA) +
  amp_cols +
  new_scale_colour() +
  geom_sf(data = wamp, aes(colour = Zone), alpha = 4/5, fill = NA) +
  wamp_col +
  plotlims +
  thetheme +
  labs(colour = NULL)

p3 <- ggplot() +
  geom_raster(data = pred_df, aes(x, y, fill = Archetype3)) +
  scale_fill_viridis() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.2) +
  geom_sf(data = amp_npz, aes(colour = ZoneName), alpha = 4/5, fill = NA) +
  amp_cols +
  new_scale_colour() +
  geom_sf(data = wamp, aes(colour = Zone), alpha = 4/5, fill = NA) +
  wamp_col +
  plotlims +
  thetheme +
  labs(colour = NULL)

p4 <- ggplot() +
  geom_raster(data = pred_df, aes(x, y, fill = Archetype4)) +
  scale_fill_viridis() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.2) +
  geom_sf(data = amp_npz, aes(colour = ZoneName), alpha = 4/5, fill = NA) +
  amp_cols +
  new_scale_colour() +
  geom_sf(data = wamp, aes(colour = Zone), alpha = 4/5, fill = NA) +
  wamp_col +
  plotlims +
  thetheme +
  labs(colour = NULL)

# make combination plot - this will take ages.
(p1 + p2) / 
  (p3 + p4)

# export - very slow!
ggsave("figures/archetype_predictions.png", 
       width = 10, height = 6, dpi = 300)
