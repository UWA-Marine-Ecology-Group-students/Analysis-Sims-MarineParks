
# Hayley's beautiful predict script - 7 covariates / 5 archetypes - enjoy!

library(raster)
library(sp)

w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
r.dir <- paste(w.dir, "SAMs/rasters", sep='/')

A_model_7_5 <- readRDS(paste(w.dir, "A_model_7_5.rds", sep='/'))
A_model_7_5

#covariates <- readRDS(paste(w.dir, "covariates.rds", sep='/'))

# load predictors --
b <- raster(paste(r.dir, "SwC_Multibeam.tiff", sep='/'))
d <- stack(paste(r.dir, "Multibeam_derivatives.tif", sep='/'))
n <- read.csv(paste(r.dir, "names.bathy.ders.csv", sep='/'))
n
names(d) <- n[,2]

d2 <- stack(d$depth * -1, d$slope, d$aspect, d$flowdir, d$tpi)   

# temp --
t1 <- raster(paste(r.dir, "SSTmean_SSTARRS.tif", sep='/'))
t2 <- raster(paste(r.dir, "SSTtrend_SSTARRS.tif", sep='/'))

t <- stack(t1, t2)
t2 <- disaggregate(t, 7.924524) #not sure what this number is? Hayley 
t3 <- resample(t2, d2)


# stack preds --
preds <- stack(d2, t3)
names(preds)

covariates <- as.data.frame(preds, xy = TRUE) 
#dim(covariates)
#head(covariates)
names(covariates) <- c('x', 'y', 'depth', 'slope', 'aspect', 'flowdir', 'tpi', 'SSTmean', 'SSTtrend')
str(covariates)
any(is.na(covariates$slope))
length(which(is.na(covariates$slope)))
covariates <- na.omit(covariates)


# predict ##
ptest2 <- predict(
  object=A_model_7_5,
  #sp.boot,
  #nboot = 100,
  newdata=covariates[,c(3:9)],
  #alpha = 0.95,
  #mc.cores = 3,
  prediction.type = "archetype"
)

#plot(ptest2)
#class(ptest2)
#str(ptest2)
#head(ptest2)

Apreds <- ptest2
head(Apreds)

SAMpreds <- cbind(covariates, Apreds)
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

plot(rasterFromXYZ(cbind(covariates[,1:2],ptest2)))
plot(A1preds)
plot(A2preds)
plot(A3preds)
plot(A4preds)
plot(A5preds)

#Response models 
par(mfrow=c(2,2))
eff.df <- effectPlotData(focal.predictors = c("depth","slope", "aspect", "tpi", "flowdir", "SSTtrend", "SSTmean"), mod = A_model_7_5)
plot(x = eff.df, object = A_model_7_5)
