# ===============================================================================
#
# Las matrix transformation
#
# ===============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 16 Feb 2020
# Last commit: 29 July 2020
# 
# Status: Finalized
#
# This file created as part of 2019 Pepperwood UAS study.
#
# ===============================================================================
#
# Description:
#
# Performs a spatial transformation on a las point cloud based on given transformation 
# matrix. The main utility of this function is to apply a transformation matrix to
# register one las data source to another. This can be computed using external
# software such as Lidar360. Function takes in a las file name (not loaded into memory) 
# and a 4x4 transformation matrix as specified below. 
#
# Prior to transformtion, converts point cloud to arbitrary coordinate system 
# whereby the point cloud's minimum (X, Y, Z) are converted to (0,0,0). Only been 
# tested on UTM projections. Needs testing on lat, long projections before use.
#
# ===============================================================================
# 
# User inputs:
#
# las_file = .las or .laz point cloud file name (.las for faster operation)
# t_matrix = 4x4 transformation matrix (see specifications below)
#
# ===============================================================================
#
# Transformation matrix:
#
# A 4x4 matrix of constants. Transformation is achieved via the following operation:
# 
# [a1 a2 a3 a4    [x
#  b1 b2 b3 b4  x  y
#  c1 c2 c3 c4     z
#  0  0  0  1]     1]
#
# such that x' = [x*a1 + y*a2 + z*a3 + 1*a4] etc. Values on the diagonal are often 
# one. Values a1 through c3 form rotational matrix while a4 to b4 shift point cloud
# in space.
#
# The following can ease creation of the transformation matrix, t. The values given
# return the same coordinates as input if transformation is performed.
#
# t <- matrix(c(
#     1.000,  0.000,  0.000,  0.000,
#     0.000,  1.000,  0.000,  0.000,
#     0.000,  0.000,  1.000,  0.000,
#     0.000,  0.000,  0.000,  1.000
#   ), nrow = 4, byrow = TRUE)

# ===============================================================================
# 
# Package dependences: 
#
# lidR, tidyverse
# 
# ===============================================================================
#
# Known problems:
#
# none
#
# ===============================================================================

suppressPackageStartupMessages(library('lidR'))
suppressPackageStartupMessages(library('tidyverse'))

lastransformation <- function(las_file, t_matrix) {
  
# =================== Convert to arbitrary coordinate system ====================
  
  crd <- readLAS(las_file, select = '')
  
  pzero = c(bbox(crd)[,1], min(crd$Z))
  
  crd <- as.matrix(crd@data) %>%
    sweep(2, pzero, '-') %>%
    t() %>%
    rbind(rep(1, ncol(.)))
  
# =========================== Perform transformation ============================
  
  crd <- t_matrix%*%crd
  
# =================== Convert back to input coordinate system ===================
  
  crd <- crd[1:3,] %>%
    t() %>%
    sweep(2, pzero, '+')
  
# ======================= Update las with new coordinates =======================
  
  las <- readLAS(las_file, select = 'c0RGB')
  
  las$X <- crd[,1]
  las$Y <- crd[,2]
  las$Z <- crd[,3]
  
  return(las)
  
}
