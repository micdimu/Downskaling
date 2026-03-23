#### function 1 predict GBM model in the plots location #####

pred_point_fun <- function(year_sel, mod, ex_plot, ntre){
  
  new_Data_loc <- ex_plot |>
    as.data.frame() |>
    add_column(yearss = as.factor(as.character(year_sel))) 
  
  pred_gbm <- predict(mod, newdata = new_Data_loc, n.trees = ntre, type = "response")
  
  df_exp <- cbind.data.frame(id = ex_plot$ID, year = year_sel, predvar = pred_gbm)
  
  return(df_exp)
}

#### function 2 run GBM models year by year and predict in the plots location #####

gbm_each <- function(year_sel, dt, ex_plot){
  
  model_file <- paste0("MODEL/model_all_", year_sel, ".RDS")
  
  dt_z <- dt |>
    select(-cell) |>
    filter(yearss == year_sel)
  
  if (file.exists(model_file)) {
    message("Loading existing model: ", model_file)
    m1 <- readRDS(model_file)
  } else {
    message("Fitting new model for year ", year_sel)
    
    m1 <- gbm(
      MAT ~ .,
      data = as.data.frame(dt_z),
      distribution = "gaussian",
      n.trees = 8000,
      interaction.depth = 3,
      shrinkage = 0.01,
      bag.fraction = 0.6,
      n.minobsinnode = 10,
      cv.folds = 5,
      keep.data = FALSE,
      verbose = FALSE
    )
    
    saveRDS(m1, model_file)
    message("Model saved to: ", model_file)
  }
  
  best1 <- gbm.perf(m1, method = "cv", plot.it = FALSE)
  
  new_Data_loc <- ex_plot |>
    as.data.frame() |>
    add_column(yearss = as.factor(as.character(year_sel)))
  
  pred_gbm <- predict(m1, newdata = new_Data_loc, n.trees = best1)
  
  df_exp <- cbind.data.frame(
    id = ex_plot$ID,
    year = year_sel,
    predvar = pred_gbm
  )
  
  return(df_exp)
}

#### function 3 predict GBM model across the study areas #####

pred_st_fun <- function(year_sel, mod, cas_p, ntre){
  
  new_Data_loc <- cas_p |>
    as.data.frame(cells = T) |>
    add_column(yearss = as.factor(as.character(year_sel))) 
  
  pred_gbm <- predict(mod, newdata = new_Data_loc, n.trees = ntre)
  
  df_exp <- cbind.data.frame(id = new_Data_loc$cell, year = year_sel, predvar = pred_gbm)
  
  return(df_exp)
}

#### function 4 predict GBM models years by years across study areas #####

pred_each_fun <- function(year_sel, cas_p, ntre){
  ls <- list.files("MODEL/", full.names = T)
  
  mod <- ls[grep(pattern = as.character(1990), ls)] |> 
    readRDS() 
  
  new_Data_loc <- cas_p |>
    as.data.frame(cells = T) |>
    add_column(yearss = as.factor(as.character(year_sel))) 
  
  pred_gbm <- predict(mod, newdata = new_Data_loc, n.trees = ntre)
  
  df_exp <- cbind.data.frame(id = new_Data_loc$cell, year = year_sel, predvar = pred_gbm)
  
  return(df_exp)
}

#### function 5 combining function 3 & 4 #####

cook_studyarea <- function(studyarea, pred, prec){
  
  cas <- vect(studyarea) |> 
    terra::project(pred) |> 
    terra::buffer(1000) |> 
    terra::ext()
  
  cas_p <- pred |> 
    terra::crop(cas)
  
  r <- rast(cas_p)  
  
  #### mod all predict ####
  
  m <- readRDS("MODEL/model_all.RDS")
  
  year_sel <- seq(1970, 2022, by = 1)
  
  pred_modall <- map(year_sel, pred_st_fun, mod = m, cas_p = cas_p, ntre = 8000) 
  
  pred_modall_w <- pred_modall |>
    reduce(bind_rows) |> 
    pivot_wider(names_from = year, values_from = predvar, names_prefix = "MAT_all_Y") |> 
    select(-id)
  
  
  MAT_all <- rep(r[[1]], ncol(pred_modall_w))
  values(MAT_all) <- pred_modall_w
  names(MAT_all) <- colnames(pred_modall_w)
  
  #### mod each predict ####
  
  pred_modeach <- map(year_sel, pred_st_fun, mod = m, cas_p = cas_p, ntre = 8000)
  
  pred_modeach_w <- pred_modeach |>
    reduce(bind_rows) |> 
    pivot_wider(names_from = year, values_from = predvar, names_prefix = "MAT_each_Y") |> 
    select(-id)
  
  
  MAT_each <- rep(r[[1]], ncol(pred_modeach_w))
  values(MAT_each) <- pred_modeach_w
  names(MAT_each) <- colnames(pred_modeach_w)
  
  #### prec interpolation ####
  
  prec_down <- prec |> 
    crop(MAT_each) |> 
    resample(MAT_each, method = "lanczos")
  
  #### export ####
  
  new_rast <- c(cas_p, MAT_all, MAT_each, prec_down)
  
  return(new_rast)
  
}

#### function 6 run GBM and predict for MAT biolimatic variables #####
# current and future scenarios across the three study areas #

gbm_bioclim <- function(x, pred, stareas){
  
  Fbio1 <- rast(x) 
  
  pred_r <- pred |> 
    resample(Fbio1, method = "bilinear") 
  
  stopifnot(compareGeom(Fbio1, pred_r, stopOnError = FALSE))
  
  X <- as.data.table(as.data.frame(pred_r, cells = T)) |> 
    drop_na()
  
  Y <- as.data.table(as.data.frame(Fbio1,  cells = T)) |> 
    drop_na()
  
  dt <- merge(Y, X, by="cell", all.x=TRUE) 
  
  dt_y <- dt |> 
    select(-cell) |> 
    select(MAT = 1, everything())
  
  m1 <- gbm(
    MAT ~ .,
    data = as.data.frame(dt_y),
    distribution = "gaussian",
    n.trees = 8000,
    interaction.depth = 3,
    shrinkage = 0.01,
    bag.fraction = 0.6,
    n.minobsinnode = 10,
    cv.folds = 5,
    keep.data = FALSE,
    verbose = FALSE
  )
  
  Fnames <- word(x, -2, sep = fixed(".tif")) |> 
    word(-1, sep = fixed("/"))
  
  saveRDS(m1, paste("MODEL/model_bioclim_", Fnames, ".RDS", sep = ""))
  
  best1 <- gbm.perf(m1, method = "cv", plot.it = FALSE)
  
  map(stareas, function(studyarea, pred, mod, ntre){
    
    sta <- vect(studyarea) |> 
      terra::project(pred) |> 
      terra::buffer(1000) |> 
      terra::ext()
    
    pred_sta <- pred |> 
      terra::crop(sta)
    
    r <- rast(pred_sta)  
    
    new_Data_loc <- pred_sta |>
      as.data.frame(cells = T) 
    
    pred_gbm <- predict(mod, newdata = new_Data_loc, n.trees = ntre)
    
    Fdown <- rep(r[[1]], 1)
    values(Fdown) <- pred_gbm
    names(Fdown) <- Fnames
    
    studynames <- studyarea |> 
      word(-2, sep = fixed(".gpkg")) |> 
      word(-1, sep = fixed("/"))
    
    writeRaster(Fdown, paste("BIOCLIM/downsc_bio_", Fnames, "_",  studynames, ".tif", sep = ""), overwrite = TRUE)
    
  }, pred = pred, mod = m1, ntre = best1)
}

#### function 7 run lanzos interpolatin for AP biolimatic variables #####
# current and future scenarios across the three study areas #

interp_bioclim <- function(x, pred, stareas){
  
  Fbio12 <- rast(x)
  
  Fnames <- word(x, -2, sep = fixed(".tif")) |> 
    word(-1, sep = fixed("/"))
  
  Fbio12_R <- Fbio12 |> 
    resample(pred, method = "lanczos")
  
  map(stareas, function(studyarea, rast_l, Fnames){
    
    sta <- vect(studyarea) |> 
      terra::project(rast_l) |> 
      terra::buffer(1000) |> 
      terra::ext()
    
    Fbio12_C <- rast_l |> 
      crop(sta)
    
    studynames <- studyarea |> 
      word(-2, sep = fixed(".gpkg")) |> 
      word(-1, sep = fixed("/"))
    
    writeRaster(Fbio12_C, paste("BIOCLIM/downsc_bio_", Fnames, "_",  studynames, ".tif", sep = ""), overwrite = TRUE)
    
  }, rast_l = Fbio12_R, Fnames = Fnames)
  
}
