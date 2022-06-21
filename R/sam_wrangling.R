###
# Project: Hayley ecomix paper
# Data:    BRUVS
# Task:    Clean wrangling script
# author:  Kingsley Griffin, Hayley Sims
# date:    May-June 2022
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
names(mb_terrain)[1] <- "mb_depth"
plot(mb_terrain)
rm(detre, zstar)

# save raster layers
writeRaster(mb_terrain, 'output/derived.tif', bylayer = TRUE, suffix = 'names')

bruvdat <- as.data.frame(raster::extract(mb_terrain, bruvdat, sp = T))          # extract terrain covariates at bruv sites
head(bruvdat)

bruvdat <- bruvdat %>% 
  mutate_at(vars(sample, scientific, family, genus, species, dataset, 
                 location, status, site), list(as.factor))
summary(bruvdat)

# clean rows with NA's in covariates
bruvdat <- bruvdat[!is.na(bruvdat$mb_depth) == TRUE, ] %>%
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

bruvdat <- bruvdat[bruvdat$scientific %in% c(sp_th), ]                          # include only sp above threshold in df
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
colnames(cov_m) <- c("sample" , "mb_depth",  "tpi", "slope", "detrend")     
summary(cov_m)        
cov_m <- cov_m[duplicated(cov_m) == FALSE, ]                                    # collapse to one row per sample
cov_m[, 2:5] <- sapply(cov_m[, 2:5], as.numeric)     
head(cov_m)

# join and setup for speciesmix
allmat   <- cbind(spw_m, cov_m)    
saveRDS(allmat, 'output/species_covariate_ecomix.rds')
