##########################################
### Import NSIDC sea ice concentration ###
##########################################

# Download sea ice concentration data from NSIDC ftp site
# Data are available either monthly or daily from 1971 to today

# Harp seal data spans 1995 to 2018
# To simplify analysis split harp seals into 4 seasons and 5 year bins

# Sea ice concentration needs to be aggregated into 20 raster layers (i.e. 4 x 5)

# Load in monthly data for each 5 year bin and output seasonal averages

# Define seasons e.g.
# Dec 1994 - Feb 1995
# Mar 1995 - May 1995
# Jun 1995 - Aug 1995
# Sep 1995 - Nov 1995

# Load libraries
require(raster)
require(tidyverse)
require(sf)

# Write the process as a function...
monthly_ice_average <- function(starts, ends){
  
  # Use RCurl library to query FTP server and list files
  url = "ftp://anonymous:wjg5@sidads.colorado.edu/DATASETS/NOAA/G02135/north/monthly/geotiff/"
  
  # Time series is split into monthly folders
  dates <- c("12_Dec", "01_Jan", "02_Feb", "03_Mar", "04_Apr", "05_May", "06_Jun",
             "07_Jul", "08_Aug", "09_Sep", "10_Oct", "11_Nov")
  
  # Generate empty raster stack to save average ice data in
  av_ice <- raster::stack()
  
  # Define progress bar
  pb1 <- txtProgressBar(min = 1, max = length(dates), style = 3) 
  
  for (j in 1:length(dates)){
    
    setTxtProgressBar(pb1, j) # update progress bar
    
    fn <- RCurl::getURL(paste0(url, dates[j], "/"),
                        ftp.use.epsv = FALSE,
                        dirlistonly = TRUE,
                        verbose = F)
  
    fn <- paste(paste0(url, dates[j], "/"),
                strsplit(fn, "\r*\n")[[1]],
                sep = "")
    
    # Only access files for sea ice concentration
    fn = fn[grep("concentration", fn)]
  
    # Find years of interest
    years <- as.numeric(substr(fn, 94, 97))
  
    # If December use start-1
    if (j == 1){
      fn <- fn[years >= (starts-1) & years <= (ends-1)]
      }
    else {
      fn <- fn[years >= starts & years <= ends]
    }
    
    # Loop through years and load ice data into stack
    x <- raster::stack()
    
    for (i in 1:length(fn)){
      ice <- raster::raster(fn[i])
      x <- raster::stack(x, ice)
    }
    raster::projection(x) = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"
    
    # 0 is ocean; 2510 pole hole; 2530 coast line; 2540 land; 2550 missing
    # 0-1000 so divide by 10 to get percentage
    x[x>1000] <- NA
    x <- x/10
    
    # Average across raster stack
    x <- raster::mean(x, na.rm = T)
    names(x) <- substr(dates[j], 4,6)
    av_ice <- raster::stack(av_ice, x)
  }
  
  # Average seasons within raster stack based on indices. See? stackApply
  indices <- c(1,1,1,2,2,2,3,3,3,4,4,4)
  av_ice <- raster::stackApply(av_ice, indices, fun = mean, na.rm = T)
  names(av_ice) <- c("Dec_Feb", "Mar_May", "Jun_Aug", "Sep_Nov")
  
  # Output the raster stack
  return(av_ice)
}


# Download data for each time period
av_ice_95 <- monthly_ice_average(starts = 1995, ends = 1999)
av_ice_00 <- monthly_ice_average(starts = 2000, ends = 2004)
av_ice_05 <- monthly_ice_average(starts = 2005, ends = 2009)
av_ice_10 <- monthly_ice_average(starts = 2010, ends = 2014)
av_ice_15 <- monthly_ice_average(starts = 2015, ends = 2019)

# Save output
raster::writeRaster(av_ice_95, filename = "Seasonal_NSIDC_Dec94_Nov99", bandorder='BIL', overwrite = T)
raster::writeRaster(av_ice_00, filename = "Seasonal_NSIDC_Dec99_Nov04", bandorder='BIL', overwrite = T)
raster::writeRaster(av_ice_05, filename = "Seasonal_NSIDC_Dec04_Nov09", bandorder='BIL', overwrite = T)
raster::writeRaster(av_ice_10, filename = "Seasonal_NSIDC_Dec09_Nov14", bandorder='BIL', overwrite = T)
raster::writeRaster(av_ice_15, filename = "Seasonal_NSIDC_Dec14_Nov19", bandorder='BIL', overwrite = T)

# ends
