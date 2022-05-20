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

# bring in and tidy bruv data
bruvdat <- read.csv("data/2020_south-west_stereo-BRUVs.complete.maxn.csv") %>%
  mutate_at(vars(sample, scientific, family, genus,  
                 species, status, site, dataset), list(as.factor)) %>% 
  glimpse()
coordinates(bruvdat) <- ~longitude + latitude 

# bring in multibeam, crop to project area and reproject to lat long
mb <- raster("data/rasters/SwC_Multibeam.tiff")
projext <- extent(288664.7 , 311265.2 , 6220416 , 6234275 )
mb      <- crop(mb, projext)
mb      <- projectRaster(mb, crs = latlon_crs)
plot(mb)
plot(bruvdat, add = TRUE)                                                       # check alignment with bruvs

# calculate and extract derivatives
mb_terrain <- terrain(mb, c('slope', "tpi"), directions = 8)                    # calculate slope and tpi (can add back in roughness, flow direction or aspect if you want)
plot(mb_terrain)

zstar <- st_as_stars(mb)                                                        # convert the bathy raster to a stars object
detre <- detrend(zstar, parallel = 8)                                           # calculate detrended bathy on your stars raster
detre <- as(object = detre, Class = "Raster")                                   # convert it back to regular raster
names(detre) <- c("detrended", "lineartrend")                                   # set the names of the new raster
plot(detre)

mb_terrain <- stack(mb, mb_terrain, detre[[1]])
names(mb_terrain)[1] <- "depth"
plot(mb_terrain)
rm(detre, zstar)

bruvdat <- as.data.frame(raster::extract(mb_terrain, bruvdat, sp = T))          # extract terrain covariates at bruv sites
head(bruvdat)

bruvdat <- bruvdat %>% 
  mutate_at(vars(sample, scientific, family, genus, species, dataset, 
                           location, status, site), list(as.factor))
summary(bruvdat)

# clean rows with NA's in covariates
bruvdat <- bruvdat[!is.na(bruvdat$depth.1) == TRUE, ] %>%
  droplevels()
summary(bruvdat)
head(bruvdat)

# Remove rare species
tot_bruv  <- length(levels(bruvdat$sample))
min_bruv  <- round(tot_bruv * 0.1)
min_bruv                                                                        # threshold number of samples

# sort out which species occur less/more than threshold
spw   <- dcast(bruvdat, sample ~ scientific, value.var = "maxn", 
               fun.aggregate = sum, drop = TRUE)                                # species long to wide format
spw[ , -1] <- sapply(spw[ , -1], as.numeric)                                    # make observations numeric
sp_pa <- ifelse(spw[, 2:ncol(spw)] >= 1, 1, 0)                                  # make presence-absence df (only species maxn)
sp_th <- colnames(sp_pa)[colSums(sp_pa) >= min_bruv]                            # list of species that occur more often than threshold

bruvdat <- bruvdat[bruvdat$scientific %in% c(sp_th), ]                         # include only sp above threshold in df
spw <- spw[ , colnames(spw) %in% c("sample", sp_th)]                            # same for wide species data
unique(bruvdat$scientific)                                                      # species remaining

# prepare species and covariates in matrix for ecomix
# species matrix first
spw_m <- spw
colnames(spw_m) <- c("sample", paste0("spp", 1:length(sp_th)))                  # long names don't work with ecomix
sp_names        <- as.data.frame(cbind(sp_th, paste0("spp", 1:length(sp_th))))  # save names for later
head(sp_names)

# covariate matrix
colnames(bruvdat)
cov_m <- bruvdat[ , c(1, 26:29)]
colnames(cov_m) <- c("sample" , "depth",  "tpi", "slope", "detrend")     
summary(cov_m)        
cov_m <- cov_m[duplicated(cov_m) == FALSE, ]                                    # collapse to one row per sample
cov_m[, 2:5] <- sapply(cov_m[, 2:5], as.numeric)     
head(cov_m)

# join and setup for speciesmix
allmat   <- cbind(spw_m, cov_m)     
sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:31),
                                                    collapse = ','),
                                     ") ~ poly(depth, 2) +
                                     poly(slope, 2) +
                                     poly(tpi, 2) +
                                     poly(detrend, 2)"))
sp_form <- ~1

set.seed(42)

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

# Look at the partial responses
par(mfrow = c(2, 2))
eff.df <- effectPlotData(focal.predictors = c("depth",  "slope", 
                                              "tpi", "detrend"),
                         mod = A_model)
plot(x = eff.df, object = A_model, se = TRUE)

# setup newdata for prediction, keeping values within range of the observations
modrange <- apply(allmat[, c('depth', 'slope', 'tpi', 'detrend')], 2, range)
newobs   <- as.data.frame(mb_terrain, na.rm = TRUE, xy = TRUE)
newobs   <- newobs %>% relocate("slope", .after = depth)
summary(newobs)
modrange
newcovs <- newobs[, -c(1, 2)]
idx <- sapply(1:ncol(newcovs), 
              function(x)ifelse(newcovs[, x] >= modrange[1, x] & 
                                  newcovs[, x] <= modrange[2, x], 1, NA))
id_inrange <- which(complete.cases(idx))
newobs     <- newobs[id_inrange, ]
colnames(newobs)[6] <- "detrend"
head(newobs)
rm(newcovs, idx, id_inrange, mb_terrain, mb)

newobs_a <- newobs[1:(nrow(newobs)/2), ]
newobs_b <- newobs[(1 + (nrow(newobs)/2)):nrow(newobs), ]

# predict
p_modela <- predict(
  object = A_model,
  #sp.boot,
  #nboot = 100,
  newdata = newobs_a, #d3
  #alpha = 0.95,
  mc.cores = 8,
  prediction.type = "archetype"
)

head(p_modela)

p_modelb <- predict(
  object = A_model,
  #sp.boot,
  #nboot = 100,
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
