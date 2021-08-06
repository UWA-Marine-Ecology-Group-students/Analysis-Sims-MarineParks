
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

model.dir <- "/homevol/anitasgiraldo2021/Analysis-Sims-MarineParks"
A_model_7_5 <- readRDS(paste(model.dir, "A_model_7_5.rds", sep='/'))
A_model_7_5

covariates <- readRDS(paste(model.dir, "covariates.rds", sep='/'))

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

plot(ptest2)
class(ptest2)
str(ptest2)
head(ptest2)

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
