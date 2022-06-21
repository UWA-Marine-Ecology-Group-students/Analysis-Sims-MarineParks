###
# Project: Hayley ecomix paper
# Data:    BRUVS
# Task:    Clean Modelling script
# author:  Kingsley Griffin, Hayley Sims
# date:    May 2022
##

library(raster)
library(dplyr)
library(stars)
library(starsExtra)
library(reshape2)
library(ecomix)

# please set your working directory manually to Analysis-Sims-MarineParks before you proceed with this script

# spatial setup
latlon_crs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# read in prepared data
allmat   <- readRDS('output/species_covariate_ecomix.rds')

# set model parameters
set.seed(42)
sp_form  <- ~1
sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:31),
                                                    collapse = ','),
                                     ") ~ poly(mb_depth, 2) +
                                     poly(slope, 2) +
                                     poly(tpi, 2) +
                                     poly(detrend, 2)"))

A_model <- species_mix(
  archetype_formula = sam_form,
  species_formula = sp_form, 
  all_formula = NULL,
  data = allmat,
  nArchetypes = 4,
  family = "negative.binomial",
  control = list(),
)

## Check model fit ----
BIC(A_model) 
AIC(A_model)
# print(A_model)

# Look at the partial responses, saving plot for publication
png('figures/archetype_response.png')
par(mfrow = c(2, 2))
eff.df <- effectPlotData(focal.predictors = c("mb_depth", "slope", 
                                              "tpi", "detrend"),
                         mod = A_model)
plot(x = eff.df, object = A_model, se = TRUE)
dev.off()

# # can also look at species
# plot(x = eff.df, object = A_model, se = TRUE, response.var = "Species")

# load raster layers
terrain_layers    <- list.files('output', '*derived', full.names = TRUE)
mb_terrain        <- stack(terrain_layers)
names(mb_terrain) <- gsub("derived_", "", names(mb_terrain))
# names(mb_terrain)[1] <- "depth.1"
plot(mb_terrain)

# setup newdata for prediction, keeping values within range of the observations
modrange <- apply(allmat[, c('mb_depth', 'slope', 'tpi', 'detrend')], 2, range)
newobs   <- as.data.frame(mb_terrain, na.rm = TRUE, xy = TRUE)
newobs   <- newobs %>% relocate("detrended", .after = tpi)                      # match order of columns to model
colnames(newobs)[6] <- "detrend"
colnames(newobs)
modrange <- round(modrange, 2)                                                  # includes a few more points
newcovs  <- newobs[, -c(1, 2)]
summary(newcovs)
modrange
idx <- sapply(1:ncol(newcovs), 
              function(x)ifelse(newcovs[, x] >= modrange[1, x] & 
                                  newcovs[, x] <= modrange[2, x], 1, NA))
id_inrange <- which(complete.cases(idx))
newobs     <- newobs[id_inrange, ]
head(newobs)
rm(newcovs, idx, id_inrange, mb_terrain)

# split the newdata to lighten load on computer
newobs_a <- newobs[1:(nrow(newobs)/2), ]
newobs_b <- newobs[(1 + (nrow(newobs)/2)):nrow(newobs), ]

# calculate variance matrix
# A_model$vcov <- vcov(A_model)

# predict
p_modela <- predict(
  object = A_model,
  # sp.boot,
  # nboot = 100,
  newdata = newobs_a, #d3
  #alpha = 0.95,
  mc.cores = 8,
  prediction.type = "archetype"
)

head(p_modela)

p_modelb <- predict(
  object = A_model,
  # sp.boot,
  # nboot = 100,
  newdata = newobs_b, #d3
  #alpha = 0.95,
  mc.cores = 8,
  prediction.type = "archetype"
)

cov_preds_a <- cbind(newobs_a, p_modela)
cov_preds_b <- cbind(newobs_b, p_modelb)
cov_preds   <- rbind(cov_preds_a, cov_preds_b)
head(cov_preds)

# grid the predictions into rasters for each archetype
coordinates(cov_preds) <- ~ x + y
A1 <- cov_preds[, 5]
A2 <- cov_preds[, 6]
A3 <- cov_preds[, 7]
A4 <- cov_preds[, 8]

gridded(A1) <- TRUE
gridded(A2) <- TRUE
gridded(A3) <- TRUE
gridded(A4) <- TRUE

p_arch <- stack(raster(A1), raster(A2), raster(A3), raster(A4))
rm(A1, A2, A3, A4)

plot(p_arch)

writeRaster(p_arch, "output/predicted.tif", bylayer = TRUE, 
            overwrite = TRUE, suffix = 'names')

