library(terra)

chelsa_download <- function(
    var = c("bio01", "bio12"),
    period,
    outdir = "CHELSA",
    gcm = NULL,
    ssp = NULL,
    overwrite = FALSE,
    quiet = TRUE
) {
  var <- match.arg(var)
  base <- "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim"
  
  # Validate var
  if (!var %in% c("bio01", "bio12")) {
    stop("This function supports only 'bio01' and 'bio12'.")
  }
  
  is_hist <- period == "1981-2010"
  
  # Build URL + local file name
  if (is_hist) {
    url <- sprintf("%s/%s/%s/CHELSA_%s_%s_V.2.1.tif", base, var, period, var, period)
    dest <- file.path(outdir, var, period, sprintf("CHELSA_%s_%s_V.2.1.tif", var, period))
  } else {
    if (is.null(gcm) || is.null(ssp)) {
      stop("For future periods you must provide both 'gcm' and 'ssp'.")
    }
    gcm_lower <- tolower(gcm)
    
    url <- sprintf(
      "%s/%s/%s/%s/%s/CHELSA_%s_%s_%s_%s_V.2.1.tif",
      base, var, period, gcm, ssp,
      gcm_lower, ssp, var, period
    )
    
    dest <- file.path(
      outdir, var, period, gcm, ssp,
      sprintf("CHELSA_%s_%s_%s_%s_V.2.1.tif", gcm_lower, ssp, var, period)
    )
  }
  
  # Download
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(dest) && !overwrite) {
    if (!quiet) message("File exists, skipping: ", dest)
    return(invisible(list(url = url, dest = dest, downloaded = FALSE)))
  }
  
  ok <- TRUE
  tryCatch(
    {
      download.file(url, destfile = dest, mode = "wb", quiet = quiet)
    },
    error = function(e) {
      ok <<- FALSE
      warning("Download failed: ", url, "\n", conditionMessage(e))
    }
  )
  
  invisible(list(url = url, dest = dest, downloaded = ok))
}


# Storico bio01
chelsa_download(var="bio01", period="1981-2010")

# Futuro bio01
chelsa_download(
  var="bio01",
  period="2011-2040",
  gcm="GFDL-ESM4",
  ssp="ssp370"
)

# Futuro bio12
chelsa_download(
  var="bio12",
  period="2011-2040",
  gcm="GFDL-ESM4",
  ssp="ssp370"
)

vars    <- c("bio01", "bio12")
periods_future <- c("2011-2040", "2041-2070")
ssps    <- c("ssp126", "ssp245", "ssp370", "ssp585")
gcms <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0") 

grid_hist <- expand_grid(var = vars, period = "1981-2010") |>
  mutate(gcm = NA_character_, ssp = NA_character_)

grid_fut <- expand_grid(
  var = vars,
  period = periods_future,
  gcm = gcms,
  ssp = ssps
)

grid <- bind_rows(grid_hist, grid_fut)

res <- pmap_dfr(grid,
  ~ chelsa_download(var = ..1, period = ..2, gcm = ..3, ssp = ..4,
                  outdir = "CHELSA", overwrite = T,
                  quiet = TRUE)
)

#### crop and overwrite ###

coords_buff <- vect("input/bbox_buffer50_lon_lat.gpkg") |> 
  terra::ext()

list.files("CHELSA", pattern = "\\.tif$", recursive = TRUE, full.names = TRUE) |>
  walk(\(f, e)
       writeRaster(
         crop(rast(f), e),
         f,
         overwrite = TRUE
       ),
       e = coords_buff)
