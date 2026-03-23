library(tidyverse)
library(terra)
library(data.table)
library(gbm)


pred <- rast("prin_downsca/pred_stack.tif")

bio_ap <- list.files("CHELSA/bio12/", full.names = T, recursive = T)
bio_mat <- list.files("CHELSA/bio01/", full.names = T, recursive = T)
stdudy_areas <- list.files("input/spatial_boundary/", full.names = T) 

map(bio_mat, gbm_bioclim, pred = pred, stareas = stdudy_areas)
map(bio_ap, interp_bioclim, pred = pred, stareas = stdudy_areas)

