# ==============================================================================
#
# LAS Ground Classification
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 10 Aug 2021
# Last commit: 10 Aug 2021
#
# Status: Completed
#
# Originated with 2019 Pepperwood UAS study. Finalized for 2021 UAS biomass study.
#
# ==============================================================================
#
# Description:
#
# Classifies ground points within LAS point clouds using CSF algorithm and the
# set of optimized parameters identified during the 2019 Pepperwood UAS study.
# First clips the uas las file to the flip box boundary.
# 
# Exports the full classified las file and a las file containing only the ground
# points. 
#
# Optimized for batch processing of large number of files using glue and foreach
#
# ==============================================================================
#
# User inputs:
#
# las_folder = Folder containing raw las files from pix4d
# spec_folder = Folder containing spectral files from pix4d. Names must be able
#       to match to las files by campaign and zone numbers
# shp_folder = Folder containing shapefiles for each uas zone. Names must be able
#       to match to las files by campaign and have zone number attribute.
# full_export = Folder for export of full las
# grnd_export = Folder for export of las containing only ground points
# grndpt_csv_export = File name for export of csv containing number of ground 
#       points identified in each UAS flight zone
#
# ==============================================================================
#
# Package dependencies:
#
# tidyverse, ggplot2, glue, 
#
# ==============================================================================
#
# Known problems:
#
# ==============================================================================

library(lidR)
library(tidyverse)
library(glue)
library(sf)
library(doParallel)

# ================================= User inputs ================================

las_folder <- 'data/las/uas/raw'

spec_folder <- 'data/spectral'

shp_gdb <- 'data/boundaries/ssu_3dforests.gdb'
shp_layer <- 'uas_zones'

full_export <- 'data/las/uas'
grnd_export <- 'data/las/icp_registration/uas_grnd_las'

grndpt_csv_export <- 'data/las/icp_registration/n_grndpts_temp.csv'



# ==============================================================================
# ========================== Classify ground points ============================
# ==============================================================================

las_files <- list.files(las_folder, pattern = 'raw')
las_file <- las_files[65]

# -------------------------- Setup cluster processing --------------------------


cl <- makeCluster(1)
registerDoParallel(cl)

grnd_pts <- foreach (
  las_file = las_files,
  .combine = 'rbind',
  .packages = c('lidR', 'tidyverse', 'glue', 'sf'),
  .export = c('las_folder', 'full_export', 'grnd_export')
) %dopar% {
  
  # -------------------------- Read in matching data ---------------------------
  
  las <- readLAS(glue('{las_folder}/{las_file}'), select = '')
  
  campaign <- str_extract(las_file, '(?<=_c)[:digit:]')
  zone <- str_extract(las_file, '(?<=_z)[:digit:]+')
  
  red <- list.files(
    spec_folder,
    pattern = glue('c{campaign}_z{zone}_red.tif'),
    full.names = TRUE
  ) %>%
    raster()
  
  green <- list.files(
    spec_folder,
    pattern = glue('c{campaign}_z{zone}_green.tif'),
    full.names = TRUE
  ) %>%
    raster()
  
  blue <- list.files(
    spec_folder,
    pattern = glue('c{campaign}_z{zone}_blue.tif'),
    full.names = TRUE
  ) %>%
    raster()
  
  rededge <- list.files(
    spec_folder,
    pattern = glue('c{campaign}_z{zone}_rededge.tif'),
    full.names = TRUE
  ) %>%
    raster()
  
  nir <- list.files(
    spec_folder,
    pattern = glue('c{campaign}_z{zone}_nir.tif'),
    full.names = TRUE
  ) %>%
    raster()
  
  shp_file <- st_read(shp_gdb, shp_layer) %>%
    filter(zone == !!zone,
           campaign == !!campaign) %>%
    st_transform(crs(las)) %>%
    st_zm() # drop Z value from polygon, produces error in clipping
  
  # ------------------------- Merge spectral data ------------------------------
  
  las <- las %>%
    filter_duplicates()
  
  las <- las %>%
    merge_spatial(source = blue,
                  attribute = 'blue') %>%
    add_lasattribute(name = 'blue', desc = 'blue')
  
  las <- las %>%
    merge_spatial(source = green,
                  attribute = 'green') %>%
    add_lasattribute(name = 'green', desc = 'green')
  
  las <- las %>%
    merge_spatial(source = red,
                  attribute = 'red') %>%
    add_lasattribute(name = 'red', desc = 'red')
  
  las <- las %>%
    merge_spatial(source = rededge,
                  attribute = 're') %>%
    add_lasattribute(name = 're', desc = 're')
  
  las <- las %>%
    merge_spatial(source = nir,
                  attribute = 'nir') %>%
    add_lasattribute(name = 'nir', desc = 'nir')
  
  ndvi <- (nir - red) / (nir + red)
  
  las <- las %>%
    merge_spatial(source = ndvi,
                  attribute = 'ndvi') %>%
    add_lasattribute(name = 'ndvi', desc = 'ndvi')
  
  ndre <- (nir - rededge) / (nir + rededge)
  
  las <- las %>%
    merge_spatial(source = ndre,
                  attribute = 'ndre') %>%
    add_lasattribute(name = 'ndre', desc = 'ndre')
  
  gndvi <- (nir - green) / (nir + green)
  
  las <- las %>%
    merge_spatial(source = gndvi,
                  attribute = 'gndvi') %>%
    add_lasattribute(name = 'gndvi', desc = 'gndvi')
  
  
  # -------------------------- Ground classification ---------------------------
  
  las <- classify_ground(
    las = las,
    algorithm = csf(
      class_threshold = 0.01,
      cloth_resolution = 0.45,
      rigidness = 3,
      time_step = 0.58,
      iterations = 500L,
      sloop_smooth = FALSE
    ),
    last_returns = FALSE
  )
  
  las@data <- las@data %>%
    mutate(Classification = replace(Classification, ndvi > 0.55, 1L))
  
  las <- las %>%
    clip_roi(shp_file)
  
  writeLAS(las,
           glue("{full_export}/{str_replace(las_file, 'raw', 'clsfd')}"))
  
  las <- filter_ground(las)
  
  if (nrow(las@data) != 0) {
    writeLAS(las,
             glue("{grnd_export}/{str_replace(las_file, 'raw', 'grnd')}"))
  }
  
  # ----------- Create data frame recording number of ground points ------------
  
  df <- tibble(campaign = campaign,
               zone = zone,
               n_ground = nrow(las@data))
  
}

write_csv(grnd_pts, grndpt_csv_export)

stopCluster(cl)

# ==============================================================================