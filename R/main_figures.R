###
# Project: Hayley ecomix paper
# Data:    BRUVS
# Task:    Make publication ready plots
# author:  Kingsley Griffin, Hayley Sims
# date:    June 2022
##

library(raster)
library(reshape2)
library(ggplot2)
library(viridis)
library(patchwork)

# get predicted archetype rasters
predictions  <- list.files("output", "*predicted", full.names = TRUE)
pred_rasters <- stack(predictions)
plot(pred_rasters)

# get NPZ boundary

# prepare predictions in long data frame for ggplot
pred_df            <- as.data.frame(pred_rasters, xy = TRUE, na.rm = TRUE)
colnames(pred_df)  <- gsub("predicted_", "", colnames(pred_df))                 # removing 'predicted'
pred_long          <- melt(pred_df, id.vars = c(1, 2))
head(pred_df)
head(pred_long)

# develop panes for each archetype
# this is slow because rasters are large and detailed
thetheme <- theme_minimal() + 
  theme(axis.title = element_blank())

p1 <- ggplot(pred_df, aes(x, y, fill = Archetype1)) +
  geom_raster() +
  scale_fill_viridis() +
  coord_equal() +
  thetheme
p1

p2 <- ggplot(pred_df, aes(x, y, fill = Archetype2)) +
  geom_raster() +
  scale_fill_viridis() +
  coord_equal() +
  thetheme
p2

p3 <- ggplot(pred_df, aes(x, y, fill = Archetype3)) +
  geom_raster() +
  scale_fill_viridis() +
  coord_equal() +
  thetheme
p3

p4 <- ggplot(pred_df, aes(x, y, fill = Archetype4)) +
  geom_raster() +
  scale_fill_viridis() +
  coord_equal() +
  thetheme
p4

# make combination plot
(p1 + p2) / 
  (p3 + p4)

# export - very slow!
ggsave("figures/archetype_predictions.png", 
       width = 10, height = 6, dpi = 300)
