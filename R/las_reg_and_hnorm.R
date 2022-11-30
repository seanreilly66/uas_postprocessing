# ==============================================================================
#
# UAS LAS Registration and Height Normalization
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 10 Aug 2021
# Last commit: 30 Nov 2022
#
# Status: Finalzied
#
# Originated with 2019 Pepperwood UAS study. Finalized for 2021 UAS biomass study.
#
# ==============================================================================
#
# Description:
#
# Performs registration of UAS LAS file using a transformation matrix produced
# by lidar360. Height normalizes the resulting transformed matrix using the
# previously generated ALS DTM.
#
# ==============================================================================
#
# User inputs:
#
# uas_las_folder = Folder containing classified uas las files. See naming 
#                  convention documentation for additional details
# als_dtm_folder = Folder containing ALS DTM raster files 
# icp_matrix_folder = Folder containing icp matrix excel files. See naming
#                     convention documentation for additional details
# grndpts_output = Output csv filename for quality control ground point check
# las_transform_path = File path to las_transformation.R function script
# cluster_size = Number of files to process at once
# 
#
# ==============================================================================
# 
# File naming convention:
#
# 1. Campaign and zone numbers
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
# 2. LAS Processing phase
#
# This script outputs records processing steps by replacing a portion of the
# filename with the completed processing step. This step of the processing
# requires the input files to contain "clsfd" in the name as output by
# uas_ground_classification.R in order to correctly identify which las files
# are in the appropriate input state and format
# 
# This is replaced by "reg_{type}" for the registered file output and hnorm
# for the height normalized final output. {Type} refers to canopy or ground 
# registration. See below for additional details
#
# 3. ICP matrices
# 
# ICP matrices must be stored in an excel file by campaign with the campaign
# number in the filename following the c## convention as described above. The 
# filename must also contain either "grnd" or "cnpy" depending on whether the
# matrices were generated using ground points only (grnd) or if canopy points
# also had to be included in instances with insufficient ground points (cnpy).
# This information is extracted from the title and appended to the output file
# name. For sites that require a mix of methods, two separate excel files must
# be generated - one "grnd" and one "cnpy" containing the subset of zones on 
# which the designated method was applied. Script will fail if a zone is 
# contained in both excel files.
# 
# Within the excel file, each zone must be stored on a separate sheet named by 
# the zone number, in the format z##. The 4x4 matrix values must be contained in 
# the first rows and columns without any titles.
# 
# ==============================================================================
#
# Package dependencies:
#
# tidyverse, glue, lidR, terra, RCSF, doParallel
#
# ==============================================================================
#
# 
# ==============================================================================
#
# Package dependencies:
#
# ==============================================================================
#
# Known problems:
#
# ==============================================================================

library(lidR)
library(tidyverse)
library(glue)
library(terra)
library(doParallel)
library(readxl)

# ================================= User inputs ================================

uas_las_folder <- 'data/las/uas_processed'

als_dtm_folder <- 'data/dtm'

icp_matrix_folder <- 'data/icp_registration/icp_matrices'

grndpts_output <- glue('data/icp_registration/registered_grndpts_{format(Sys.time(), "%Y%m%d_%H%M")}.csv')

las_transform_path <- 'R/las_transformation.R'

cluster_size <- 10

# ==============================================================================
# =================== Registration and height normalization ==================== 
# ==============================================================================

# -------------- Setup parallel processing and load matching files ------------- 

clsfd_files <- list.files(path = uas_las_folder, pattern = 'clsfd')

cl <- makeCluster(cluster_size)
registerDoParallel(cl)

grndpts <- foreach (
  uas_file = clsfd_files,
  .combine = 'rbind',
  .packages = c('lidR', 'tidyverse', 'glue', 'terra', 'readxl')
) %dopar% {
  
  source(las_transform_path)
  
  campaign <- str_extract(uas_file, '(?<=c)[:digit:]+')
  zone <- str_extract(uas_file, '(?<=z)[:digit:]+')
  
  als_dtm <- list.files(als_dtm_folder,
                        pattern = glue('c{campaign}.+\\.tif'),
                        full.names = TRUE) %>%
    str_subset(glue('z{zone}')) %>%
    rast()
  
  icp_matrix <- list.files(path = icp_matrix_folder,
                           pattern = glue('c{campaign}'),
                           full.names = TRUE)
  
  # Identify matrix for sites that required canopy alignment at some plots
  if (length(icp_matrix > 1)) {
    for (x in icp_matrix) {
      sheets = excel_sheets(x)
      if (glue('z{zone}') %in% sheets) {
        icp_matrix <- x
      }
    }
  }
  
  if (str_detect(icp_matrix, 'grnd')) {
    matrix_type <- 'grnd'
  } else if (str_detect(icp_matrix, 'cnpy')) {
    matrix_type <- 'cnpy'
  } else {
    stop('Incorrect icp matrix file designation: ', icp_matrix)
  }

  icp_matrix <- read_xlsx(
    path = icp_matrix,
    sheet = glue('z{zone}'),
    col_names = FALSE,
    range = 'A1:D4'
  ) %>%
    as.matrix()

  # ------------------------- ICP matrix registration --------------------------

  las <- lastransformation(las_file = glue(uas_las_folder, '/', uas_file),
                           t_matrix = icp_matrix)
  
  reg_filename <- str_replace(uas_file, 'clsfd', glue('reg_{matrix_type}')) %>%
    glue(uas_las_folder, '/', .)

  writeLAS(las, file = reg_filename)

  # ---------- Extract ground points and dtm Z for offset confirmation ---------

  grnd <- las %>%
    filter_ground() %>%
    merge_spatial(source = als_dtm,
                  attribute = 'dtm_z')

  # --------------------- DTM based height normalization ---------------------

  las <- normalize_height(las, als_dtm, na.rm = TRUE, add_lasattribute = TRUE)

  hnorm_filename <- str_replace(uas_file, 'clsfd', 'hnorm') %>%
    glue(uas_las_folder, '/', .)
  
  writeLAS(las, file = hnorm_filename)

  # ---------------------- Output ground points dataset ----------------------

  if (matrix_type == 'grnd') {
    grnd <- grnd@data %>%
      rename(uas_z = Z) %>%
      select(uas_z, dtm_z) %>%
      add_column(zone = zone,
                 campaign = campaign)
  } else {
    grnd <- data.frame(
      uas_z = NA,
      dtm_z = NA,
      zone = zone,
      campaign = campaign
    )
  }
  
}

write_csv(grndpts, grndpts_output)

stopCluster(cl)