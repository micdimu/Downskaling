library(tidyverse)
library(terra)
library(sf)
library(data.table)
library(gbm)

#### Load coordinates and buffer ####
coords <- list.files(path = "input/", pattern = "plots", full.names = T) |>
  map(vect) |> 
  (\(.) do.call(rbind, .))() 

coords_buff <- vect("input/bbox_buffer50_lon_lat.gpkg") |> 
  terra::project("EPSG:32632")

#### DEM crop ####

dem <- rast("input/demIta_10m_crop.tif") |> 
  terra::crop(coords_buff) |> 
  aggregate(
  fact = 10,
  fun  = mean,
  na.rm = TRUE
)

##### 1 Slope & Aspect (Eastness & Northness) #####

slope_deg <- terrain(
  dem,
  v        = "slope",   
  unit     = "degrees", 
  neighbors = 8
)

names(slope_deg) <- "slope_deg"

aspect_deg <- terrain(
  dem,
  v        = "aspect",
  unit     = "degrees",
  neighbors = 8
)

names(aspect_deg) <- "aspect_deg"

# eastness and northness

asp <- terrain(dem, "aspect", unit = "radians")
eastness  <- sin(asp)
northness <- cos(asp)

names(eastness) <- "eastness"
names(northness) <- "northness"



##### 2 TRI (Terrain Ruggedness Index) #####

tri <- terrain(dem, v = "TRI")
names(tri) <- "TRI"

rough <- terrain(dem, v = "roughness")
names(rough) <- "roughness"


##### 3 TPI (Topographic Position Index) #####

tpi <- terrain(dem, v = "TPI")
names(tpi) <- "TPI" 


##### 4 SD elevation #####

w3 <- matrix(1, nrow = 3, ncol = 3)

sd_elev_3x3 <- focal(
  dem,
  w   = w3,
  fun = sd,
  na.rm = TRUE
)

names(sd_elev_3x3) <- "sd_elev_3x3"

##### 5 Flow direction #####

flowdir <- terrain(dem, v = "flowdir")
names(flowdir) <- "flowdir"


##### 6 TWI (Topographic wetness index) #####

slope <- terrain(dem, "slope", unit = "radians")

w_area <- cellSize(dem, unit = "m")
flow_acc <- flowAccumulation(flowdir, weight = w_area)

flow_acc[flow_acc <= 0] <- NA
slope[slope <= 0]       <- NA

twi <- log(flow_acc / tan(slope))


#### Stack and export ####

topo_stack <- c(
  dem,          # elevazione
  slope_deg,    # pendenza (°)
  aspect_deg,   # esposizione (°)
  tri,          # Terrain Ruggedness Index
  rough,        # roughness (max-min locale)
  tpi,          # Topographic Position Index
  flowdir,      # direzione di flusso D8
  sd_elev_3x3,  # micro-rilievo (sd quota 3x3)
  twi,          # Topographic Wetness Index,
  eastness,     # Eastness
  northness     # Northness
) |> 
  project("epsg:4326")


#### Latitudinal Predictor ####

lat <- init(topo_stack, "y")

#### Solar Radiation ####

solrad <- rast("input/GHI.tif") |> 
  crop(topo_stack) |> 
  resample(topo_stack,  method = "bilinear") # GHI – Global Horizontal Irradiation

########## export ########

idy <- names(topo_stack) %in% c("dem_10m", "eastness", "northness", "slope_deg")

pred <- c(topo_stack[[idy]],
          solrad,
          lat)      # layers: DEM, eastness, northness, slope, solRad, latitude

writeRaster(pred, "input/pred_stack.tif", overwrite = TRUE)

pred_plots <- cbind.data.frame(id = coords$Plot_name, extract(pred, coords))
 
write.csv(pred_plots, "input/pred_plots.csv")

