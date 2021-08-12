## Script to extract bathymetry covariates from raster files ####
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(sp)
library(sf)
library(raster)
library(rgdal)
library(spatstat)
library(sp)
devtools::install_github('skiptoniam/ecomix')
library(ecomix)
library(reshape2)
library(tidyr)
library(devtools)
library(BBmisc) # to normalize
library(corrplot)
library(tibble)
library(rasterVis)
library(RColorBrewer)
library(broman)

# clear workspace ----
rm(list = ls())

# set working directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "H:/Github/GB_2015_Survey"
# Set data directory - to read the data from
dt.dir <- (paste(w.dir, "Data/Tidy", sep='/'))
h.dir <- (paste(w.dir, "Data/Habitat/BRUV Style annotation/tidy data"))
s.dir <- (paste(w.dir, "shapefiles", sep='/'))
# Set graph directory - to save plots
p.dir <- paste(w.dir, "Plots", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')
o.dir <- paste(w.dir, "outputs", sep='/')


# Load Multibeam ----
b <- raster(paste(r.dir, "SwC_Multibeam.tiff", sep='/'))
plot(b)

# crop to extent --
#e <- drawExtent()
e <- extent(288664.7 , 311265.2 , 6220416 , 6234275 )
b <- crop(b, e)
plot(b)
b # 4x4m resolution

#### Transform from utm to lat long ----

# open reference file 
ref <- raster(paste(r.dir, "GB-SW_250mBathy.tif", sep='/'))
ref.crs <- proj4string(ref)

b <- projectRaster(b, crs = ref.crs)
plot(b)
proj4string(b) # check it worked


# Calculate bathy derivatives ----
s <- terrain(b, 'slope')
plot(s)
a <- terrain(b, 'aspect')
plot(a)
r <- terrain(b, 'roughness')
plot(r)
t <- terrain(b, 'TPI')
plot(t)
f <- terrain(b, 'flowdir')
plot(f)

ders <- stack(b,s,a,r,t,f)
names(ders) <- c("depth", "slope",  "aspect" ,  "roughness"  ,   "tpi" ,   "flowdir")



# save stack of derivatives
writeRaster(ders, paste(r.dir, "Multibeam_derivatives.tif", sep='/'))

#Script 1 #####################################

# Load data ----

study <- "2020_south-west_stereo-BRUVs"

df <- read.csv(paste(dt.dir, paste(study, "complete.maxn.csv", sep='.'), sep = '/'))%>%
  mutate_at(vars(sample, scientific, family, genus,  species, status, site, dataset), list(as.factor)) %>% # make these columns as factors
  glimpse()
head(df)
str(df)

dfs <- df
coordinates(dfs) <- ~longitude+latitude 

# Get bathy and detrended derivatives ----
#b <- raster(paste(r.dir, "Multibeam_derivatives.tif", sep='/')) #HS did not run for new bathy

# plot bathy and points
plot(b) 
plot(dfs, pch = 20, cex = 1, add=T) 

b[b>-35] <-NA
plot(b)

# Load  derivatives ----
bds <- stack(paste(r.dir, "Multibeam_derivatives.tif", sep='/')) #changed to "Multibeam_derivatives.tif"
names(bds)
names2 <- read.csv(paste(r.dir, "names.bathy.ders.csv", sep='/'))
names(bds) <- names2$x
names(bds)
plot(bds)


# Extract bathy derivatives from data points --
dfs <- raster::extract(bds, dfs, sp = T)
str(dfs)
head(dfs)



# save maxn with covariates ----
write.csv(dfs, paste(dt.dir, "2020_sw_maxn.env-cov.csv", sep='/'))


###       ###       ###       ###


# Get SST covariates ----
t1 <- raster(paste(r.dir, "SSTmean_SSTARRS.tif", sep='/'))
t2 <- raster(paste(r.dir, "SSTsterr_SSTARRS.tif", sep='/'))
t3 <- raster(paste(r.dir, "SSTtrend_SSTARRS.tif", sep='/'))

ts <- stack(t1, t2, t3)
plot(ts)
plot(ts$SSTmean_SSTARRS)

dfs <- raster::extract(ts, dfs, sp=T)
head(dfs)


#### save mxn with bathy ders and temp ----
write.csv(dfs, paste(dt.dir, "2020_sw_maxn.env-cov.csv", sep='/'))

##### Script 3 Archetypes



# 1. Load data ----
df <- read.csv(paste(dt.dir, "2020_sw_maxn.env-cov.csv", sep = '/'))%>%
  mutate_at(vars(sample, scientific, family, genus, species, dataset, location, status, site), list(as.factor)) %>% # make these columns as factors
  # At some point filter for successful count
  glimpse()
head(df)
summary(df)
levels(df$scientific) # 140 species

## Add log depth ----
df <- df[,-36]
df$logdepth <- log((df$depth.1 * -1))
any(is.na(df$depth.1))
any(is.na(df$logdepth))
length(which(is.na(df$logdepth)))

# clean rows with NA's in covariates
summary(df)
df <- df[!is.na(df$depth.1) == TRUE, ] %>%
  droplevels()
summary(df)
head(df)


# 2. Remove sp that are encountered in less than a defined threshold percentage of samples ----
# as per Foster et al 2015 ----
tot_bruv  <- length(levels(df$sample))
min_bruv  <- round(tot_bruv * 0.025)
min_bruv                                                                        # threshold number of samples

# sort out which species occur less/more than threshold
spw   <- dcast(df, sample ~ scientific, value.var = "maxn", 
               fun.aggregate = sum, drop = TRUE)                                # species long to wide format
spw[ , -1] <- sapply(spw[ , -1], as.numeric)                                    # make observations numeric
sp_pa <- ifelse(spw[, 2:ncol(spw)] >= 1, 1, 0)                                  # make presence-absence df (only species maxn)
sp_th <- colnames(sp_pa)[colSums(sp_pa) >= min_bruv]                            # list of species that occur more often than threshold

df  <- df[df$scientific %in% c(sp_th), ]                                        # include only sp above threshold in df
spw <- spw[ , colnames(spw) %in% c("sample", sp_th)]                            # same for wide species data
unique(df$scientific)                                                           # species remaining

# 3. Prepare species and covariates in matrix for ecomix ----
spw_m <- spw
colnames(spw_m) <- c("sample", paste0("spp", 1:length(sp_th)))                  # long names don't work with ecomix
sp.names.no     <- as.data.frame(cbind(sp_th, paste0("spp", 1:length(sp_th))))  # save names for later
head(sp.names.no)

cov_m <- df[ , c(2, 27:36)]
colnames(cov_m) <- c("sample"  ,"depth",  "slope"  ,  "aspect", "roughness",
                     "tpi", "flowdir", "SSTmean", "SSTster", 
                     "SSTtrend", "logdepth")     
summary(cov_m)                                                                  # remove sstster, as it's all one value
cov_m         <- cov_m[ , -9]
#cov_m$depth   <- cov_m$depth * -1                                               # ecomix doesn't seem to like negatives?
cov_m$flowdir <- as.numeric(cov_m$flowdir)
cov_m <- cov_m[duplicated(cov_m) == FALSE, ]                                    # collapse to one row per sample
cov_m[, 2:9] <- sapply(cov_m[, 2:9], as.numeric)                                # make observations numeric

# cov_ms <- scale(cov_m[, -1])                                                    # scale/standardise covariates (leave for now)
allmat <- cbind(spw_m, cov_m)                                                   # join species and scaled covariates
str(allmat)                                                                     # checking form of each variable
# allmat[ , c(2:28, 30:38)] <- sapply(allmat[, c(2:28, 30:38)], as.numeric)       # ecomix seems to need numeric



# 4. Calculate correlation coefficient between covariates ----

## Check correlation among covariates ---



# 7. Final model ----

sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:71),
                                                    collapse = ','),
                                     ") ~ poly(depth, 2) + poly(aspect, 2) + poly(tpi, 2) +
                                     poly(flowdir, 2) + poly(SSTtrend, 2) + poly(SSTmean, 2) +
                                     poly(slope, 2)"))


sp_form <- ~1

A_model <- species_mix(
  archetype_formula = sam_form,
  #poly(slope, degree = 2, raw = TRUE) + poly(tpi, degree = 2, raw = TRUE) + poly(aspect, degree = 2, raw = TRUE) +
  #poly(temp_mean, degree = 2, raw = TRUE) + poly(temp_trend, degree = 2, raw = TRUE)),
  species_formula = sp_form, #stats::as.formula(~1),
  all_formula = NULL,
  data=allmat,
  nArchetypes = 5,
  family = "negative.binomial",
  #offset = NULL,
  #weights = NULL,
  #bb_weights = NULL,
  #size = NULL, # for presence absence - benthic point data
  #power = NULL, # for tweedie : eg. biomass data
  control = list(), # for tuning the model if needed
  #inits = NULL, # if you have fitted the model previously: use the same values
  #standardise = FALSE, # has been removed in new update it scales
  #titbits = TRUE # could turn this off
)


## Check model fit ----
BIC(A_model) # this gives a value of BIC
print(A_model)
A_model$coefs
A_model$alpha
A_model$beta
A_model$lofl
A_model$gamma
A_model$delta
A_model$theta

# look at the partial response of each covariate using:
par(mfrow=c(2,2))
eff.df <- effectPlotData(focal.predictors = c("depth","aspect", "tpi", "flowdir", "SSTtrend", "SSTmean", "slope"), mod = A_model)
plot(x = eff.df, object = A_model)


# 11. Probability of each sp. belonging to each archetype ----
arch_prob <- as.data.frame(A_model$tau)
head(arch_prob)
names(arch_prob)
head(sp.names.no)

arch <- cbind(arch_prob, sp.names.no)
head(arch)
arches <- cbind(arch_prob, sp.names.no)
head(arches)

# Get archetype with maximum probability for each sp --
for(i in 1:nrow(arches)){
  arches$arch[i] <- colnames(arches)[which.max(arches[i, 1:ncol(arch_prob)])]   # takes column name with maximum probability
  arches$p[i]    <- max(arches[i, 1:ncol(arch_prob)])                           # takes maximum probability value
}
head(arches)
ggplot(data = arch, aes(p_arch)) + geom_bar()                                   # check all archetypes have >1 member


# Probability of each sp. belonging to each archetype ----
arch_prob <- as.data.frame(A_model$tau)
head(arch_prob)
names(arch_prob)
head(sp.names.no)

arch <- cbind(arch_prob, sp.names.no)
head(arch)

# Get archetype with maximum probability for each sp --
arch2 <- arch_prob %>%
  tibble::rownames_to_column() %>% # 1. add row ids as a column
  gather(column, value, -rowname) %>%
  dplyr::group_by(rowname) %>%
  dplyr::filter(rank(-value) == 1) %>%
  glimpse()

head(arch2)
str(arch2)
arch2 <- as.data.frame(arch2)
names(arch2) <- c('new.names', 'archetype', 'prob')

colnames(sp.names.no)[2] <- c("new.names")

arch3 <- arch2 %>%
  dplyr::left_join(sp.names.no) %>%
  dplyr::mutate_at(vars(archetype), list(as.factor)) %>%
  glimpse()

str(arch3)
summary(arch3)


# to save this list --
#write.csv(arches, paste(o.dir, "species_archetypes.csv", sep ='/'))



# 12. Plots ----

plot(A_model)

preds <- c("slope","depth", "aspect", "tpi", "flowdir", "SSTtrend", "SSTmean")


ef.plot <- effectPlotData(preds, A_model)
head(ef.plot)

#ef.plot$bathy
#plot(ef.plot, A_model)


#sp.boot <- species_mix.bootstrap(
#A_model,
#nboot = 100,
#type = "BayesBoot",
#mc.cores = 2,
#quiet = FALSE
#)

#saveRDS(sp.boot, paste(model.dir, "sp.boot.rds", sep ='/'))


#plot(ef.plot,
#A_model, 
#sp.boot,
#)


# still not sure what to use this for --
#sp.var <- vcov(
# A_model,
#sp.boot = sp.boot,
#method = "BayesBoot",
#nboot = 10,
#mc.cores = 2
#)



## 13. Predict ----
#model.dir <- "/homevol/anitasgiraldo2021/Analysis-Sims-MarineParks"
#A_model <- readRDS(paste(model.dir, "A_model.rds", sep='/'))
#A_model

# load predictors --

# bathy --
b <- raster(paste(r.dir, "SwC_Multibeam.tiff", sep='/'))
plot(b)
# cut to 150 m depth --
#b[b < -150] <- NA
#b[b > -30] <- NA
#plot(b)

# derivatives --
d <- stack(paste(r.dir, "Multibeam_derivatives.tif", sep='/'))
plot(d)
names(d)
n <- read.csv(paste(r.dir, "names.bathy.ders.csv", sep='/'))
n
#n$covs <- c("depth", "slope", "aspect", "roughness", "tpi", "flowdir")
names(d) <- n[,2]

# crop bathy to stack --
#b2 <- crop(b, d)
#d <- mask(d, b2)
#plot(d)

# stack preds --
#d2 <- stack(d$depth, d$slope)
#plot(d2)
d2 <- stack(d$depth, d$slope, d$aspect, d$flowdir, d$tpi)    

# make depth positive (matching the conversion we did before modelling)
#plot(d2)
#plot(d2$depth)

#e <- drawExtent()
###e <- extent(114.7212, 114.9377, -34.13335, -34.12439)

###d2.2 <- crop(d2, e)
###plot(d2.2)

###d3 <- as.data.frame(d2.2, xy = TRUE)

# temp --
t1 <- raster(paste(r.dir, "SSTmean_SSTARRS.tif", sep='/'))
t2 <- raster(paste(r.dir, "SSTtrend_SSTARRS.tif", sep='/'))

t <- stack(t1, t2)

t2 <- disaggregate(t, 7.924524)
t3 <- resample(t2, d2)

#plot(t3)
#plot(t2)


# crop bathy to stack --
#t2 <- crop(t, b)
#d <- mask(d, b2)
#plot(d)

# stack preds --
preds <- stack(d2, t3)
names(preds)
plot(preds)
plot(d)


d3 <- as.data.frame(preds, xy = TRUE)
dim(d3)
head(d3)
names(d3) <- c('x', 'y', 'depth', 'slope', 'aspect', 'flowdir', 'tpi', 'SSTmean', 'SSTtrend')
str(d3)

any(is.na(d3$slope))
length(which(is.na(d3$slope)))
d3 <- na.omit(d3)
str(d3)


# predict ##
ptest2 <- predict(
  object=A_model,
  #sp.boot,
  #nboot = 100,
  newdata=d3[,c(3:9)],
  #alpha = 0.95,
  #mc.cores = 3,
  prediction.type = "archetype"
)


plot(rasterFromXYZ(cbind(d3[,1:2],ptest2)))

Apreds <- ptest2
head(Apreds)

SAMpreds <- cbind(d3, Apreds)
head(SAMpreds)

coordinates(SAMpreds) <- ~x+y
A1 <- SAMpreds[,8]
A2 <- SAMpreds[,9]
A3 <- SAMpreds[,10]
A4 <- SAMpreds[,11]
A5 <- SAMpreds[,12]

gridded(A1) <- TRUE
gridded(A2) <- TRUE
gridded(A3) <- TRUE
gridded(A4) <- TRUE
gridded(A5) <- TRUE

A1preds <- raster(A1)
A2preds <- raster(A2)
A3preds <- raster(A3)
A4preds <- raster(A4)
A5preds <- raster(A5)

Allpreds <- stack(A1preds, A2preds, A3preds, A4preds, A5preds)
plot(Allpreds)

plot(rasterFromXYZ(cbind(d3[,1:2],ptest2)))
plot(A1preds)
plot(A2preds)
plot(A3preds)
plot(A4preds)
plot(A5preds)

#Response models again
par(mfrow=c(2,2))
eff.df <- effectPlotData(focal.predictors = c("depth","slope", "aspect", "tpi", "flowdir", "SSTtrend", "SSTmean"), mod = A_model)
plot(x = eff.df, object = A_model)




##########################################################
