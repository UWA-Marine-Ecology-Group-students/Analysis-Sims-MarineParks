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



## so far just using bathy covariates ##

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
min_bruv  <- round(tot_bruv * 0.05)
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
cov_m$depth   <- cov_m$depth * -1                                               # ecomix doesn't seem to like negatives?
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

sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:49),
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
  nArchetypes = 6,
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
eff.df <- effectPlotData(focal.predictors = c("depth","slope", "aspect",  "flowdir", "SSTtrend", "SSTmean"), mod = A_model)
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

# to save this list --
#write.csv(arches, paste(o.dir, "species_archetypes.csv", sep ='/'))



# 12. Plots ----

plot(A_model)

preds <- c("slope","depth", "aspect", "flowdir", "SSTtrend", "SSTmean")


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
d2 <- stack(d$depth * -1, d$slope, d$aspect, d$flowdir)   

# make depth positive (matching the conversion we did before modelling)
plot(d2)

plot(d2$depth)

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

plot(t3)
plot(t2)


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
names(d3) <- c('x', 'y', 'depth', 'slope', 'aspect', 'flowdir', 'SSTmean', 'SSTtrend')
str(d3)
any(is.na(d3$slope))
length(which(is.na(d3$slope)))
d3 <- na.omit(d3)
str(d3)

any(is.na(d3$SSTmean))
length(which(is.na(d3$SSTmean)))
d3 <- na.omit(d3)
str(d3)



#memory.size(30000)

# predict ##
ptest2 <- predict(
  object=A_model,
  #sp.boot,
  #nboot = 100,
  newdata=d3[,c(3:8)],
  #alpha = 0.95,
  #mc.cores = 3,
  prediction.type = "archetype"
)

plot(ptest2)
plot(rasterFromXYZ(cbind(d3[,1:2],ptest2)))
class(ptest2)
str(ptest2)
head(ptest2)
#ptest2$ptPreds

Apreds <- ptest2
head(Apreds)

#d3.1<- d3[ , c(2, 27:36)]

SAMpreds <- cbind(d3, Apreds)
head(SAMpreds)


coordinates(SAMpreds) <- ~x+y
A1 <- SAMpreds[,3]
A2 <- SAMpreds[,4]
A3 <- SAMpreds[,5]
A4 <- SAMpreds[,6]

gridded(A1) <- TRUE
gridded(A2) <- TRUE
gridded(A3) <- TRUE
gridded(A4) <- TRUE

A1preds <- raster(A1)
A2preds <- raster(A2)
A3preds <- raster(A3)
A4preds <- raster(A4)

Allpreds <- stack(A1preds, A2preds, A3preds, A4preds)
plot(Allpreds)


##########################################################


#Temp 

ders <- stack(b,s,a,r,t,f)
names(ders) <- c("depth", "slope",  "aspect" ,  "roughness"  , "flowdir")




# PREDICT using test model ----


#sp.boot <- species_mix.bootstrap(
 # test_model,
  #nboot = 100,
  #type = "BayesBoot",
  mc.cores = 10,
  quiet = FALSE
)




# load predictors --

# bathy --
b <- stack(paste(r.dir, "SW_bathy.derivatives-to-260m.tif", sep='/'))
plot(b)
names(b)
nam <- read.csv(paste(r.dir, "names.bathy.ders.csv", sep='/'))
names(b) <- nam$x

d<- dropLayer(b, c(2:6))
d
plot(d)
e <- drawExtent()
d <- crop(d, e)

# cut to 150 m depth --
values(d)[values(d) < -150] = NA
values(d)[values(d) > -30] = NA
#b[b < -150] <- NA
#b[b > -30] <- NA

s <- crop(b$slope, d)
s <- mask(s, d)

a <- crop(b$aspect, d)
a <- mask(a, d)

t  <- crop(b$tpi, d)
t <- mask(t, d)

f <- crop(b$flowdir, d)
f <- mask(f, d)

ds <- stack(d,s,a,t,f)
plot(ds)


# temp --
t1 <- raster(paste(r.dir, "SSTsterr_SSTARRS.tif", sep='/'))
t2 <- raster(paste(r.dir, "SSTtrend_SSTARRS.tif", sep='/'))

t <- stack(t1, t2)

t2 <- disaggregate(t, 7.924524)
t3 <- resample(t2, d)

plot(t3)


# crop bathy to stack --
#t2 <- crop(t, b)
#d <- mask(d, b2)
#plot(d)

# stack preds --
preds <- stack(d, t3)
names(preds)
plot(preds)
plot(d)


pr <- as.data.frame(preds, xy = TRUE)
dim(pr)
head(pr)
names(pr) <- c('x', 'y', 'depth', 'slope', 'aspect', 'tpi',  'flowdir', 'SSTster', 'SSTtrend')
str(pr)
any(is.na(pr$slope))
length(which(is.na(pr$slope)))
pr <- na.omit(pr)
str(pr)
head(pr)


# predict ##
ptest2 <- predict(
  test_model,
  sp.boot,
  #nboot = 100,
  pr[,c(3:9)],
  #alpha = 0.95,
  mc.cores = 10,
  prediction.type = "archetype"
)


class(ptest2)
str(ptest2$ptPreds)
head(ptest2$ptPreds)
ptest2$ptPreds

Apreds <- ptest2$ptPreds
head(Apreds)

SAMpreds <- cbind(pr, Apreds)
head(SAMpreds)


coordinates(SAMpreds) <- ~x+y
A1 <- SAMpreds[,8]
A2 <- SAMpreds[,9]
A3 <- SAMpreds[,10]
#A4 <- SAMpreds[,7]

gridded(A1) <- TRUE
gridded(A2) <- TRUE
gridded(A3) <- TRUE
#gridded(A4) <- TRUE

A1preds <- raster(A1)
A2preds <- raster(A2)
A3preds <- raster(A3)
#A4preds <- raster(A4)

Allpreds <- stack(A1preds, A2preds, A3preds)
plot(Allpreds)

writeRaster(Allpreds, paste(o.dir, "pred-3a-7cov.tif", sep='/'))



# PREDICT using test model b ----


sp.boot <- species_mix.bootstrap(
  test_model_b,
  nboot = 100,
  type = "BayesBoot",
  mc.cores = 10,
  quiet = FALSE
)

# vairance -- Still not sure how to do this an the summary
sp.var <- vcov(
  test_model_b,
  object2 = NULL,
  method = "BayesBoot",
  nboot = 10,
  mc.cores = 10
)


summary.species_mix(test_model_b, sp.var)

# load predictors --

# bathy --
b <- stack(paste(r.dir, "SW_bathy.derivatives-to-260m.tif", sep='/'))
plot(b)
names(b)
nam <- read.csv(paste(r.dir, "names.bathy.ders.csv", sep='/'))
names(b) <- nam$x

d<- dropLayer(b, c(2:6))
d
plot(d)
e <- extent(114.4438, 115.0793, -34.28812, -33.58538)
d <- crop(d, e)
plot(d)

# cut to 150 m depth --
values(d)[values(d) < -143] = NA
values(d)[values(d) > -35] = NA
#b[b < -150] <- NA
#b[b > -30] <- NA

s <- crop(b$slope, d)
s <- mask(s, d)

a <- crop(b$aspect, d)
a <- mask(a, d)

t  <- crop(b$tpi, d)
t <- mask(t, d)

f <- crop(b$flowdir, d)
f <- mask(f, d)

ds <- stack(d,a,f,s)
plot(ds)


# temp --
t1 <- raster(paste(r.dir, "SSTsterr_SSTARRS.tif", sep='/'))
t2 <- raster(paste(r.dir, "SSTtrend_SSTARRS.tif", sep='/'))

tx <- stack(t1, t2)
plot(tx)

ta <- disaggregate(t1, 7.924524)
tb <- disaggregate(t2, 7.924524)
ta <- resample(ta, ds)
plot(ta)

tb <- resample(tb, ds)
plot(tb)

t <- stack(ta, tb)


# crop bathy to stack --
#t2 <- crop(t, b)
#d <- mask(d, b2)
#plot(d)

# stack preds --
preds <- stack(ds, t)
names(preds)
plot(preds)


pr <- as.data.frame(preds, xy = TRUE)
dim(pr)
head(pr)
names(pr) <- c('x', 'y', 'depth','aspect', 'slope', 'flowdir', 'SSTster', 'SSTtrend')
str(pr)
any(is.na(pr$slope))
length(which(is.na(pr$slope)))
pr <- na.omit(pr)
str(pr)
head(pr)


# predict ##
ptest2 <- predict(
  test_model_b,
  sp.boot,
  #nboot = 100,
  pr[,c(3:8)],
  #alpha = 0.95,
  mc.cores = 10,
  prediction.type = "archetype"
)


class(ptest2)
str(ptest2$ptPreds)
head(ptest2$ptPreds)
ptest2$ptPreds

Apreds <- ptest2$ptPreds
head(Apreds)

SAMpreds <- cbind(pr, Apreds)
head(SAMpreds)


coordinates(SAMpreds) <- ~x+y
head(SAMpreds)
A1 <- SAMpreds[,7]
A2 <- SAMpreds[,8]
A3 <- SAMpreds[,9]
#A4 <- SAMpreds[,7]

gridded(A1) <- TRUE
gridded(A2) <- TRUE
gridded(A3) <- TRUE
#gridded(A4) <- TRUE

A1preds <- raster(A1)
A2preds <- raster(A2)
A3preds <- raster(A3)
#A4preds <- raster(A4)

Allpreds <- stack(A1preds, A2preds, A3preds)
plot(Allpreds)

writeRaster(Allpreds, paste(o.dir, "Ecomix5.tif", sep='/'), overwrite =T)

test <- log(Allpreds)
plot(test)



# PREDICT USING STAND COVS ----

## Standarize the covariates --
d3.stand <- BBmisc::normalize(d3[,c(3:5)], method = "standardize", range = c(0,1))
colnames(d3.stand) <- colnames(d3[,c(3:5)])
dim(d3.stand)
head(d3.stand)


ptest2 <- predict(
  A_model,
  sp.boot,
  d3.stand,
  prediction.type = "archetype"
)

class(ptest2)
str(ptest2$ptPreds)
head(ptest2$ptPreds)

Bpreds <- ptest2$ptPreds
head(Bpreds)

SAMpreds <- cbind(d3, Bpreds)
head(SAMpreds)


coordinates(SAMpreds) <- ~x+y
A1 <- SAMpreds[,4]
A2 <- SAMpreds[,5]
A3 <- SAMpreds[,6]
A4 <- SAMpreds[,7]

gridded(A1) <- TRUE
gridded(A2) <- TRUE
gridded(A3) <- TRUE
gridded(A4) <- TRUE

A1preds <- raster(A1)
A2preds <- raster(A2)
A3preds <- raster(A3)
A4preds <- raster(A4)

AllpredsB <- stack(A1preds, A2preds, A3preds, A4preds)
plot(AllpredsB)

#####################
p1 <- stack(paste(o.dir, "pred-3a-5cov-notpiaspect.tif", sep='/'))
p2 <- stack(paste(o.dir, "pred-3a-6cov-notpi.tif", sep='/'))
p3 <- stack(paste(o.dir, "pred-3a-7cov.tif", sep='/'))
p4 <- stack(paste(o.dir, "pred-3a-notpiSSTtrend.tif", sep='/'))

plot(p4)
names(p4) <- c("Archetype1",  "Archetype2", "Archetype3")

# define breaks manually
breaks1 <- c(-Inf, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 1000, Inf)
p <- levelplot(Allpreds, par.settings=RdBuTheme(), at=breaks1, 
               colorkey=list(height=0.8, labels=list(at=breaks1, labels=round(breaks1, 2))))



p


## FOR MAPS ####

ab <- readOGR(paste(s.dir, "Australiaboundary67.shp", sep='/'))

amp <- readOGR(paste(s.dir, "AustraliaNetworkMarineParks.shp", sep='/'))

wamp <- readOGR(paste(s.dir, "All_StateReserves.shp", sep='/'))

plot(test$Archetype3, main = "Archetype 3")
plot(ab, add=T)
plot(amp, add=T)
plot(wamp, add=T)

dev.off()

p     


## Get quantiles ----
A1 <- A1preds
A2 <- A2preds
A3 <- A3preds

A1[A1>5000] <- 5000
plot(A1)

A2[A2>5000] <- 5000
plot(A2)

A3[A3>1000] <- 1000
plot(A3)

plot(Allpreds$Archetype1)
plot(Allpreds$Archetype2)
plot(Allpreds$Archetype3)
plot(test$Archetype1)
ecomix.quant <- c(0,0.05,0.95,1)
ecomix.cuts <- quantile(Allpreds$Archetype1, ecomix.quant)
ecomix.cuts
brk<- quantile(c(Allpreds$Archetype2))
#catA1 <- cut(Allpreds$Archetype1, breaks=ecomix.cuts, na.rm=TRUE)
#plot(catA1)
plot(Allpreds$Archetype2)
ahist<-hist(Allpreds$Archetype2,
            breaks= 100,
            main="Histogram Archetype 2",
            col="blue",  # changes bin color
            xlab= "Abundance")
ahist
breaks1 <- c(0,1,10,50,100,1000,Inf)
breaks2 <- c(0, 1, 10, 50, 100, 1000, 5000, 38405.51)
breaks2 <- c(0, 1, 10, 50, 100, 1000, 2500, 5000)
p <- levelplot(Allpreds$Archetype2, par.settings=RdBuTheme(), at=breaks1, 
               colorkey=list(height=0.8, labels=list(at=breaks1, labels=round(breaks1, 2))))



p


plot_crayons()
yel <- brocolors("crayons")["Canary"]
or <- brocolors("crayons")["Vivid Tangerine"]
g0 <- brocolors("crayons")["Electric Lime"]
g1 <- brocolors("crayons")["Screamin Green"]
g2 <- brocolors("crayons")["Sea Green"]
b1 <- brocolors("crayons")["Caribbean Green"]
b2 <- brocolors("crayons")["Blue Green"]
b3 <- brocolors("crayons")["Navy Blue"]
b4 <- brocolors("crayons")["Violet Blue"]
r1 <- brocolors("crayons")["Radical Red"]
r2 <- brocolors("crayons")["Maroon"]

#pal <- colorRampPalette(c("red","blue", "green"))
pal <- colorRampPalette(c(yel, g0, g1, b1, g2, b2, b3, r1))

plot(A2,
     #Allpreds$Archetype2,
     breaks = breaks2,
     #col = terrain.colors(6),
     col = pal(7),
     main = "Archetype 2",
     legend = FALSE)
#r.range <- c(minValue(Allpreds$Archetype2), maxValue(Allpreds$Archetype2))
#plot(Allpreds$Archetype2)
plot(A2,
     #Allpreds$Archetype2,
     legend.only = TRUE,
     #col = terrain.colors(6),
     col = pal(7),
     legend.width=1, legend.shrink=0.75,
     axis.args = list(at = c(714, 1428, 2142, 2856, 3570, 4284, 5000), 
                      labels = c('1','10','50','100','1000', '2500','5000')))
#axis.args=list(at =  c(0, 1, 10, 50, 100, 1000, 5000, Inf)),
#labels = breaks2), 
#             cex.axis=0.6),
#legend.args=list(text='Abundance', side=1, font=2, line=2.5, cex=0.8))

dev.off()

plot(ab, col = 'orange', alpha= 0.5, add=T)
plot(amp, add=T)
plot(wamp, add=T)

## USING ggplot ----
r = Allpreds$Archetype2 #raster object
#preparing raster object to plot with geom_tile in ggplot2
r_points = rasterToPoints(r)
r_df = data.frame(r_points)
head(r_df) #breaks will be set to column "layer"
r_df$cuts=cut(r_df$Archetype2, breaks=c(1, 10, 50, 100, 1000, 5000, 39000)) #set breaks

ggplot(data=r_df) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  #scale_fill_brewer("Legend_title",type = "seq", palette = "Greys") +
  scale_fill_manual("Abundance",  values = pal(7)) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude")






