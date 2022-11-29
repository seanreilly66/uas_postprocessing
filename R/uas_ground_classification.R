# ==============================================================================
#
# UAS Ground Classification and Spectral Merge
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 10 Aug 2021
# Last commit: 28 Nov 2022
#
# Status: Completed
#
# Originated with 2019 Pepperwood UAS study. Finalized for 2021 UAS biomass study.
#
# ==============================================================================
#
# Description:
#
# Performs initial postprocessing steps on Pix4D exports. 
# 
# Merges UAS las file and spectral data
#
# Classifies ground points within LAS point clouds using CSF algorithm and the
# set of optimized parameters identified during the 2019 Pepperwood UAS study.
# 
# Exports the full classified las file and a las file containing only the ground
# points as required for icp registration. 
#
# Optimized for batch processing of large number of files using glue and foreach
#
# ==============================================================================
#
# User inputs:
#
# raw_folder = Folder containing raw las files from pix4d
# spec_folder = Folder containing spectral files from pix4d
# shp_folder = Folder containing shapefiles for each uas zone. Names must be able
#       to match to las files by campaign and have zone number attribute.
# full_export = Folder for export of full las
# grnd_export = Folder for export of las containing only ground points
# grndpt_csv_export = File name for export of csv containing number of ground 
#       points identified in each UAS flight zone
# cluster_size = Number of files to process at once
#
# ==============================================================================
# 
# File naming convention:
# 
# This script uses a naming convention based on campaign and zone numbers to
# match files when batch processing. All files must contain the following:
# 
# campaign number in the format c## 
# zone or plot number in the format z##
# 
# example for uas file from campaign 1 zone 3: c1_z3
# 
# campaign and zone numbers do not need to be adjacent and can be separated by
# other information (e.g., c1_ebr2_z5) so long as the name does not contain
# other instances of c## or z##
#
# 
# Spectral files need to also contain the corresponding band spelled out in 
# their file name in lower case. These need to be:
# red, green, blue, rededge, nir
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

raw_folder <- 'data/las/uas_raw'

spec_folder <- 'data/spectral'

full_export <- 'data/las/uas_processed'
grnd_export <- 'data/icp_registration/uas_ground'

grndpt_csv_export <- 'data/icp_registration/n_grndpts_temp.csv'

cluster_size <- 10

# ==============================================================================
# ========================== Classify ground points ============================
# ==============================================================================

raw_files <- list.files(raw_folder)

# -------------------------- Setup cluster processing --------------------------


cl <- makeCluster(cluster_size)
registerDoParallel(cl)

grnd_pts <- foreach (
  las_file = raw_files,
  .combine = 'rbind',
  .packages = c('lidR', 'tidyverse', 'glue', 'sf'),
  .export = c('raw_folder', 'full_export', 'grnd_export')
) %dopar% {
  
  # -------------------------- Read in matching data ---------------------------
  
  las <- readLAS(glue('{raw_folder}/{las_file}'), select = '')
  
  campaign <- str_extract(las_file, '(?<=c)[:digit:]')
  zone <- str_extract(las_file, '(?<=z)[:digit:]+')
  
  spec_files <- list.files(
    spec_folder,
    pattern = glue('c{campaign}'),
    full.names = TRUE) %>% 
    str_subset(glue('z{zone}'))
  
  red <- list.files(
    spec_folder,
    pattern = glue('c{campaign}_z{zone}_red.tif'),
    full.names = TRUE
  ) %>%
    rast()
  
  green <- list.files(
    spec_folder,
    pattern = glue('c{campaign}_z{zone}_green.tif'),
    full.names = TRUE
  ) %>%
    rast()
  
  blue <- list.files(
    spec_folder,
    pattern = glue('c{campaign}_z{zone}_blue.tif'),
    full.names = TRUE
  ) %>%
    rast()
  
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