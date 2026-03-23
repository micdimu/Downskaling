library(tidyverse)
library(terra)
library(data.table)
library(gbm)
source("Source_fun.R")

pred <- rast("input/pred_stack.tif")

clim_var <- rast("input/Clim_stack_crop.tif")
idx <- grepl("Tavg", names(clim_var))

mat <- clim_var[[idx]]
prec <- clim_var[[!idx]]

names(prec) <- names(prec) |> 
  str_replace("DownscaledPrcp", "AP_interp_Y") |>
  str_replace("YearlySum_cogeo", "")  

unduetre <- list.files("spatial_boundary/", full.names = T) |> 
  map(cook_studyarea, pred = pred, prec = prec)

terra::writeRaster(unduetre[[1]], "Down_MAT_AP_Casentino.tif", overwrite = T)
terra::writeRaster(unduetre[[2]], "Down_MAT_AP_Lagorai.tif", overwrite = T)
terra::writeRaster(unduetre[[3]], "Down_MAT_AP_ParcoSirentiVelino.tif", overwrite = T)



