# uas_postprocessing
Scripts to postprocess UAS files, including spectral merging, ALS registration, and height normalization

## Workflow

1. ALS LAS files are fed through als_dtm_generation.R to generate raster and las DTMs
2. UAS outputs from Pix4D are processed through las_ground_classification.R to merge spectral data and to classify ground points.
3. Ground point outputs from ALS and UAS are registered in Lidar360 and the resulting transformation matrices are stored in an xlsx file
4. Classified UAS files are fed into uas_reg_and_hnorm.R, along with raster DTM and transformation matrices, which registers the point cloud and height normalizes it.

## Folder structure

This repository contains the folder structure as currently expected by the contained scripts. As a result, it can be used "as is" by putting the raw files in their respective places. 

## Naming convention

The scripts use a nuanced naming convention to match different data from a campaign and zone. Additional information is contained in each of the scripts.
