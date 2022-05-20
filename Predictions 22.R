# Use A_model1.rds if A_model not in local environment 
# Final model 4 Archetypes - Depth, Slope, TIP and Detrended Bathy

# Look at the partial response
par(mfrow=c(1,1))
eff.df <- effectPlotData(focal.predictors = c("depth","slope", "detrend", "tpi"), mod = A_model)
plot(x = eff.df, object = A_model, se=TRUE)


# Derivatives --
d <- stack(paste(r.dir, "Multibeam_derivatives.tif", sep='/'))
n <- read.csv(paste(r.dir, "names.bathy.ders.csv", sep='/'))
names(d) <- n[,2]

d2 <- stack(d$depth, d$slope, d$detrend, d$tpi) 
preds <-d2
names(preds)
plot(preds)


d3 <- as.data.frame(preds, xy = TRUE)
names(d3) <- c('x', 'y', 'depth', 'slope', 'tpi', 'detrend')

any(is.na(d3$slope))
length(which(is.na(d3$slope)))
d3 <- na.omit(d3)
str(d3)

## skips code - this little bit of code will drop site which sit outside of model space.
modrange <- apply(allmat[,c('depth', 'slope', 'tpi', 'detrend')],2,range)
newobs <- d3
summary(newobs)
newobs.range <- newobs[,-1:-2]
idx <- sapply(1:ncol(newobs.range),function(x)ifelse(newobs.range[,x]>=modrange[1,x]&newobs.range[,x]<=modrange[2,x],1,NA))
env.in.range.idx <- which(complete.cases(idx))
newobs2 <- newobs[env.in.range.idx,]


# predict ##
ptest2 <- predict(
  object=A_model,
  #sp.boot,
  #nboot = 100,
  newdata=newobs2[,c(3:6)], #d3
  #alpha = 0.95,
  #mc.cores = 3,
  prediction.type = "archetype"
)

Apreds <- ptest2
head(Apreds)
SAMpreds <- cbind(newobs2, Apreds)
head(SAMpreds)


coordinates(SAMpreds) <- ~x+y
A1 <- SAMpreds[,5]
A2 <- SAMpreds[,6]
A3 <- SAMpreds[,7]
A4 <- SAMpreds[,8]

gridded(A1) <- TRUE
gridded(A2) <- TRUE
gridded(A3) <- TRUE
gridded(A4) <- TRUE

A1preds <- raster(A1)
A2preds <- raster(A2)
A3preds <- raster(A3)
A4preds <- raster(A4)

plot(A1preds)
plot(A2preds)
plot(A3preds)
plot(A4preds)


#Predict Plots 

plot(A1preds, ylab="Latitude", xlab="Longitude", main="Archetype 1 (4m x 4m)" )
plot(A2preds, ylab="Latitude", xlab="Longitude", main="Archetype 2 (4m x 4m)" )
plot(A3preds, ylab="Latitude", xlab="Longitude", main="Archetype 3 (4m x 4m)" )
plot(A4preds, ylab="Latitude", xlab="Longitude", main="Archetype 4 (4m x 4m)" )


#increase area  (x by 2.5)
aggregated.r.A1 <- raster::aggregate(A1preds, fact = 12.5, fun = sum)  
plot(aggregated.r.A1, ylab="Latitude", xlab="Longitude", main="Archetype 1 (50m x 50m)") 

aggregated.r.A2 <- raster::aggregate(A2preds, fact = 12.5, fun = sum)  
plot(aggregated.r.A2, ylab="Latitude", xlab="Longitude", main="Archetype 2 (50m x 50m)") 

aggregated.r.A3 <- raster::aggregate(A3preds, fact = 12.5 , fun = sum)  
plot(aggregated.r.A3, ylab="Latitude", xlab="Longitude", main="Archetype 3 (50m x 50m)") 

aggregated.r.A4 <- raster::aggregate(A4preds, fact = 12.5 , fun = sum)  
plot(aggregated.r.A4, ylab="Latitude", xlab="Longitude", main="Archetype 4 (50m x 50m)") 


##########################################################


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
glimpse(arch3)


# Plots ----

plot(A_model)
preds <- c("slope","depth", "aspect", "tpi", "flowdir")

####################################



# clear workspace ----
rm(list = ls())

# set working directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
dt.dir <- (paste(w.dir, "Data/Tidy", sep='/'))
h.dir <- (paste(w.dir, "Data/Habitat/BRUV Style annotation/tidy data"))
s.dir <- (paste(w.dir, "shapefiles", sep='/'))
p.dir <- paste(w.dir, "Plots", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')
o.dir <- paste(w.dir, "outputs", sep='/')


# Load Multibeam ----
b <- raster(paste(r.dir, "SwC_Multibeam.tiff", sep='/'))
plot(b)

# crop to extent --
e <- extent(288664.7 , 311265.2 , 6220416 , 6234275 )
b <- crop(b, e)
plot(b)

# Transform from utm to lat long ----

# open reference file 
ref <- raster(paste(r.dir, "GB-SW_250mBathy.tif", sep='/'))
ref.crs <- proj4string(ref)
b <- projectRaster(b, crs = ref.crs)
plot(b)

# Derivatives ----
s <- terrain(b, 'slope')
a <- terrain(b, 'aspect')
r <- terrain(b, 'roughness')
t <- terrain(b, 'tpi')
f <- terrain(b, 'flowdir')

ders <- stack(b,s,a,r,t,f)
names(ders) <- c("depth", "slope",  "aspect" ,  "roughness"  ,   "tpi" ,   "flowdir")

writeRaster(ders, paste(r.dir, "Multibeam_derivatives.tif", sep='/'))

# Load data ----

study <- "2020_south-west_stereo-BRUVs"

df <- read.csv(paste(dt.dir, paste(study, "complete.maxn.csv", sep='.'), sep = '/'))%>%
  mutate_at(vars(sample, scientific, family, genus,  species, status, site, dataset), list(as.factor)) %>% 
  glimpse()

dfs <- df
coordinates(dfs) <- ~longitude+latitude 

plot(dfs, pch = 20, cex = 1, add=T) 


# Load  derivatives ----
bds <- stack(paste(r.dir, "Multibeam_derivatives.tif", sep='/')) 
names2 <- read.csv(paste(r.dir, "names.bathy.ders.csv", sep='/'))
names(bds) <- names2$x

# Extract bathy derivatives from data points --
dfs <- raster::extract(bds, dfs, sp = T)

# save maxn with covariates ----
write.csv(dfs, paste(dt.dir, "2020_sw_maxn.env-cov.csv", sep='/'))

df <- read.csv(paste(dt.dir, "2020_sw_maxn.env-cov.csv", sep = '/'))%>%
  mutate_at(vars(sample, scientific, family, genus, species, dataset, location, status, site), list(as.factor)) %>% # make these columns as factors
  # At some point filter for successful count
  glimpse()

# clean rows with NA's in covariates
summary(df)
df <- df[!is.na(df$depth.1) == TRUE, ] %>%
  droplevels()
summary(df)
head(df)


# Remove sp that are encountered in less than a defined threshold percentage of samples ----
# as per Foster et al 2015 ----
tot_bruv  <- length(levels(df$sample))
min_bruv  <- round(tot_bruv * 0.1)
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

cov_m <- df[ , c(2, 27:32)]
colnames(cov_m) <- c("sample" , "depth",  "slope"  ,  "aspect", "roughness",
                     "tpi", "flowdir", "detrend")     
summary(cov_m)        

cov_m1 <-df

cov_m         <- cov_m[ , -9]                                       
cov_m$flowdir <- as.numeric(cov_m$flowdir)
cov_m <- cov_m[duplicated(cov_m) == FALSE, ]                                    # collapse to one row per sample
cov_m[, 2:7] <- sapply(cov_m[, 2:7], as.numeric)     

#scale aspect 

cov_ms <-(cov_m[4] - min(cov_m[4]))/(max(cov_m[4] - min(cov_m[4])))


#address circularity 
x <- cov_ms
cov_msa1 <- (sin((2*pi)*(x/1))+cos((2*pi)*(x/1)))
cov_msa <-cbind(cov_m[2:3])
cov_msa2 <- cbind(cov_m[5:7])

allmat <- cbind(spw_m, cov_msa1, cov_msa, cov_msa2)     

sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:31),
                                                    collapse = ','),
                                     ") ~ poly(depth, 2) +
                                     poly(detrend, 2) +
                                     poly(tpi, 2) +
                                     poly(slope, 2)"))


sp_form <- ~1

A_model <- species_mix(
  archetype_formula = sam_form,
  species_formula = sp_form, 
  all_formula = NULL,
  data=allmat,
  nArchetypes = 4,
  family = "negative.binomial",
  control = list(),
)


## Check model fit ----
BIC(A_model) 
AIC(A_model)
print(A_model)

# Look at the partial response
par(mfrow=c(1,1))
eff.df <- effectPlotData(focal.predictors = c("depth","tpi","detrend","slope"), mod = A_model)
plot(x = eff.df, object = A_model, se=TRUE)


# Derivatives --
d <- stack(paste(r.dir, "Multibeam_derivatives.tif", sep='/'))
n <- read.csv(paste(r.dir, "names.bathy.ders.csv", sep='/'))
names(d) <- n[,2]

d2 <- stack(d$depth, d$slope, d$flowdir, d$tpi, d$aspect) 
preds <-d2
names(preds)
plot(preds)


d3 <- as.data.frame(preds, xy = TRUE)
names(d3) <- c('x', 'y', 'depth', 'slope', 'flowdir', 'tpi', 'aspect')

any(is.na(d3$slope))
length(which(is.na(d3$slope)))
d3 <- na.omit(d3)
str(d3)

#save model 
saveRDS(A_model, "A_model1.rds")
#my_model <- readRDS("A_model1.rds")

## skips code - this little bit of code will drop site which sit outside of model space.
modrange <- apply(allmat[,c('depth', 'slope', 'detrend', 'tpi')],2,range)
newobs <- d3
summary(newobs)
newobs.range <- newobs[,-1:-2]
idx <- sapply(1:ncol(newobs.range),function(x)ifelse(newobs.range[,x]>=modrange[1,x]&newobs.range[,x]<=modrange[2,x],1,NA))
env.in.range.idx <- which(complete.cases(idx))
newobs2 <- newobs[env.in.range.idx,]


# predict ##

ptest2 <- predict(
  object=A_model,
  #sp.boot,
  #nboot = 100,
  newdata=newobs2[,c(3:6)], #d3
  #alpha = 0.95,
  #mc.cores = 3,
  prediction.type = "archetype"
)

Apreds <- ptest2
head(Apreds)
SAMpreds <- cbind(d3, Apreds)
head(SAMpreds)


coordinates(SAMpreds) <- ~x+y
A1 <- SAMpreds[,5]
A2 <- SAMpreds[,6]
A3 <- SAMpreds[,7]
A4 <- SAMpreds[,8]


gridded(A1) <- TRUE
gridded(A2) <- TRUE
gridded(A3) <- TRUE
gridded(A4) <- TRUE


A1preds <- raster(A1)
A2preds <- raster(A2)
A3preds <- raster(A3)
A4preds <- raster(A4)


plot(A1preds)
plot(A2preds)
plot(A3preds)
plot(A4preds)


#Predict Plots 


plot(A1preds, ylab="Latitude", xlab="Longitude", main="A1 - 4m x 4m" )
aggregated.r.A1 <- raster::aggregate(A1preds, fact = 10, fun = sum)  
plot(aggregated.r.A1, ylab="Latitude", xlab="Longitude", main="A1") 

plot(A2preds, ylab="Latitude", xlab="Longitude", main="A2 - 4m x 4m" )
aggregated.r.A2 <- raster::aggregate(A2preds, fact = 10, fun = sum)  
plot(aggregated.r.A2, ylab="Latitude", xlab="Longitude", main="A2") 

plot(A3preds, ylab="Latitude", xlab="Longitude", main="A3 - 4m x 4m" )
aggregated.r.A3 <- raster::aggregate(A3preds, fact = 10 , fun = sum)  
plot(aggregated.r.A3, ylab="Latitude", xlab="Longitude", main="A3") 

plot(A4preds, ylab="Latitude", xlab="Longitude", main="A4 - 4m x 4m" )
aggregated.r.A4 <- raster::aggregate(A4preds, fact = 10 , fun = sum)  
plot(aggregated.r.A4, ylab="Latitude", xlab="Longitude", main="A4") 


##########################################################


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
glimpse(arch3)


# Plots ----

plot(A_model)
preds <- c("slope","depth", "aspect", "tpi", "flowdir")



aggregated.r.A1 <- raster::aggregate(A1preds, fact = 10, fun = sum)  
plot(aggregated.r.A1, ylab="Latitude", xlab="Longitude", main="A1 40m x 40m") 

aggregated.r.A2 <- raster::aggregate(A2preds, fact = 10, fun = sum)  
plot(aggregated.r.A2, ylab="Latitude", xlab="Longitude", main="A2 40m x 40m") 

aggregated.r.A3 <- raster::aggregate(A3preds, fact = 10, fun = sum)  
plot(aggregated.r.A3, ylab="Latitude", xlab="Longitude", main="A3 40m x 40m") 



aggregated.r.A8 <- raster::aggregate(A8preds, fact = 5, fun = sum)  
plot(aggregated.r.A8, ylab="Latitude", xlab="Longitude", main="A8 40m x 40m") 

aggregated.r.A2 <- raster::aggregate(A2preds, fact = 10, fun = sum)  
plot(aggregated.r.A2, ylab="Latitude", xlab="Longitude", main="A2 40m x 40m") 





