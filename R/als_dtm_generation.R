# ==============================================================================
#
# ALS DTM LAS generation
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
# Created as part of 2021-2022 UAS biomass study.
#
# ==============================================================================
#
# Description:
#
# Generates DTM from ALS las files and converts the DTM to a LAS file format for
# use in ground icp registration of UAS files.
#
# Setup for parallel batch processing of all las files in one parent folder
#
# ==============================================================================
#
# User inputs:
#
# las_folder = Folder containing las files with classified ground points
# dtm_las_output_folder = Location to which las dtms should be exported
# cluster_size = Number of files to process in parallel
# dtm_resolution = Resolution (in meters) of intermediate DTM raster
#
# ==============================================================================
#
# Package dependencies:
#
# lidr, tidyverse, glue, sf, doParallel
#
# ==============================================================================
#
# Known problems:
#
# none
#
# ==============================================================================

library(lidR)
library(tidyverse)
library(glue)
library(terra)
library(doParallel)

# ================================= User inputs ================================

las_folder <- 'data/las/als'

raster_output_folder <- 'data/dtm'

las_output_folder <- 'data/icp_registration/als_dtm_las'

cluster_size = 10

dtm_resolution = 0.5


# ==============================================================================
# ============================ Generate DTM las file =========================== 
# ==============================================================================

las_files <- list.files(las_folder, full.names = TRUE) 

cl <- makeCluster(cluster_size)
registerDoParallel(cl)

foreach (
  lf = las_files,
  .packages = c('lidR', 'tidyverse', 'glue', 'sf', 'terra')
) %dopar% {
  
  # --------------------------- Generate raster DTM ---------------------------- 
  
  dtm <- readLAS(lf, filter = '-filter_class 2') %>%
    rasterize_terrain(res = dtm_resolution, algorithm = tin())
  
  dtm_name <- lf %>%
    str_replace('\\.las$', '_dtm\\.tif') %>%
    str_replace(las_folder, raster_output_folder)
  
  writeRaster(dtm, dtm_name)

  # ----------------------------- Convert to LAS -------------------------------
  
  dtm_las <- dtm %>%
    as.data.frame(xy = TRUE) %>%
    filter(!is.na(Z),
           Z > -100) %>%
    rename(
      X = x, 
      Y = y
    ) %>%
    LAS()
  
  projection(dtm_las) <- crs(dtm)
  
  # -------------------------------- Export ------------------------------------
  
  output_name <- lf %>%
    str_replace('\\.las$', '_dtm\\.las') %>%
    str_replace(las_folder, las_output_folder)
  
  writeLAS(dtm_las, output_name)
  
}
  
stopCluster(cl)

# ==============================================================================