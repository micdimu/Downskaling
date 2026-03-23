library(tidyverse)
library(terra)
library(data.table)
library(gbm)
library(furrr)
source("Source_fun.R")
#tmpFiles(remove = TRUE)

pred <- rast("input/pred_stack.tif")

clim_var <- rast("input/Clim_stack_crop.tif")
idx <- grepl("Tavg", names(clim_var))

mat <- clim_var[[idx]]
prec <- clim_var[[!idx]]

names(prec) <- names(prec) |>
  str_replace("DownscaledPrcp", "AP_interp_Y") |>
  str_replace("YearlySum_cogeo", "")

coords <- list.files(path = "input/", pattern = "plots", full.names = T) |>
  map(vect) |> 
  (\(.) do.call(rbind, .))() 

pred_r <- pred |> 
  resample(mat, method = "bilinear") 

stopifnot(compareGeom(mat, pred_r, stopOnError = FALSE))

### function ###

years <- names(mat)   
year_num <- gsub(".*?(\\d{4}).*", "\\1", years)

X <- as.data.table(as.data.frame(pred_r, cells = T)) |> 
  drop_na()
Y <- as.data.table(as.data.frame(mat,  cells = T)) |> 
  drop_na()

N <- length(X$cell) ## cambiare N per fare un test rapido senza aspettare ore
## tempo di calcolo con tutte le celle sono 19h ##

set.seed(1234)
cells <- sample(X$cell, N)

# training X sulle celle campionate (nrow == length(cells) garantito)
X_f <- X[cell %in% cells]
Y_f <- Y[cell %in% cells]


# 2) long format: una riga = cella-anno
dt <- melt(Y_f, id.vars="cell", variable.name="yearss", value.name="MAT") |> 
  merge(X_f, by="cell", all.x=TRUE) |> 
  mutate(yearss = gsub(".*?(\\d{4}).*", "\\1", yearss)) |> 
  mutate(yearss = as.factor(yearss))

dt_y <- dt |> 
  select(-cell) 

extr_plot_loc <- terra::extract(pred, coords) |> 
  mutate(ID = coords$Plot_name)

### NOT TO RUN - VERY TIME CONSUMING ###
start <- Sys.time()

m <- gbm(
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

Sys.time() - start

#saveRDS(m, "MODEL/model_all.RDS")

m <- readRDS("MODEL/model_all.RDS")

best <- gbm.perf(m, method = "cv", plot.it = FALSE)

##### Predict ######

seq_years <- dt_y |> 
  pull(yearss) |> 
  as.character() |> 
  as.numeric() |> 
  unique()

pred_modall <- map(seq_years, pred_point_fun, mod = m, ex_plot = extr_plot_loc, ntre = best) |> 
  reduce(bind_rows) |> 
  pivot_wider(names_from = year, values_from = predvar, names_prefix = "MAT_Y")


write.csv(pred_modall, "YEARLY/MATyearly_downscaled_plots.csv")

##### Model Each and Predict ######

mod_each <- map(seq_years, gbm_each, dt = dt, ex_plot = extr_plot_loc)

mod_each_i <- mod_each |> 
  reduce(bind_rows) |> 
  pivot_wider(names_from = year, values_from = predvar, names_prefix = "MAT_Y")

write.csv(mod_each_i, file = "YEARLY/MATyearly_downscaled_plots_each.csv")

##### Annual Precipitation #####

prec_r <- prec |> 
  resample(pred, method = "lanczos")

prec_plots <- cbind.data.frame(id = coords$Plot_name, extract(prec_r, coords))

write.csv(prec_plots, "YEARLY/APyearly_downscaled_plots_each.csv")

##### check correlation ####

MATall <- read.csv("YEARLY/MATyearly_downscaled_plots.csv", row.names = 1)
MATeach <- read.csv("YEARLY/MATyearly_downscaled_plots_each.csv", row.names = 1)

MATall1 <- MATall |> 
  pivot_longer(cols = -id) |> 
  rename("all" = "value")

MATeach1 <- MATeach |> 
  pivot_longer(cols = -id) |> 
  rename("each" = "value")

checkdiff <- full_join(MATall1, MATeach1) |> 
  mutate(diff = all - each)

checkdiff |> 
  pull(diff) |> 
  hist()

cor(checkdiff$each, checkdiff$all)

modEvA::RMSE(obs = checkdiff$each, pred = checkdiff$all)

