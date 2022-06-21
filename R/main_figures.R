###
# Project: Hayley ecomix paper
# Data:    BRUVS
# Task:    Make publication ready plots
# author:  Kingsley Griffin, Hayley Sims
# date:    June 2022
##

library(raster)
library(ggplot2)

# estimated effects plot
esteff <- readRDS("output/effect_estimates.rds")

for(i in 1:4){
  idat     <- esteff[[i]]
  idf      <- as.data.frame(idat)
  idf$term <- colnames(idat)[i]
  if(i == 1){
    esteff_df <- idf
  } else {
    esteff_df <- rbind(esteff_df, idf)
  }
}

head(esteff_df)

ggplot()