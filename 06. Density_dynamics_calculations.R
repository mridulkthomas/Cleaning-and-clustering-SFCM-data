## Quantifying cell densities and biovolumes of phytoplankton communities and functional 
## groups using scanning flow cytometry, machine learning and unsupervised clustering
## Mridul K. Thomas, Simone Fontana, Marta Reyes & Francesco Pomati
#
## Code designed to be used with the dataset available at 10.5281/zenodo.977773
#
## Install all required packages before proceeding
# 
# list.of.packages <- c('data.table','randomForest', 'plot3Drgl', 'dplyr',
#                       'lubridate', 'dtplyr', 'corrplot')
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if (length(new.packages)) install.packages(new.packages)
# 
## This package has a different source and needs to be installed separately
# source("http://bioconductor.org/biocLite.R")
# biocLite("flowPeaks")


### Script 6. Calculate densities of the total phytoplankton community
# Code takes as input: 
#   1) Estimated cell densities in each measurement from the scanning flow cytometer
#   2) Cleaned data files, from Script 4
# Note: This script will therefore need to be adapted based on the formatting of the
#     FCM or SFCM instrument output
# See associated manuscript for details

# Set working directory
setwd()

# Load packages
library(data.table)
library(lubridate)

# Calculate number of particles in cleaned data files
cleanfiles <- list.files(path = './Script 4. Cleaning by clustering/output/', 
                         pattern = c('cleaned'), full.names = TRUE)

# Create empty matrix to populate with information about each cleaned file
clean <- matrix(1, nrow = length(cleanfiles), ncol = 2)
cleansummary <- as.data.table(clean)
cleansummary <- cleansummary[, V1 := as.integer(V1)]
cleansummary <- cleansummary[, V2 := as.character(V2)]
rm(clean)

# Loop through all cleaned files and record number of particles, date & time from each file

for (index in 1L:length(cleanfiles)) {
  print(index)
  cl <-       fread(cleanfiles[index], select = 1)
  set(cleansummary, i = index, j = 1L, dim(cl)[1])
  set(cleansummary, i = index, j = 2L, 
      as.character(ymd_hm(gsub('u',':', substr(cleanfiles[index], 57, 72)))))
}

setnames(cleansummary,c('count_clean','datetime_orig'))

# Add measured number of particles and particle density of raw files
raw_dens <- read.csv('./Script 6. Calculating cell density/input/raw_file_cell_density_measurements.csv') 

cleansummary$count_unclean <- raw_dens$particles_measured
cleansummary$dens_unclean <- raw_dens$density
cleansummary$filenames <- raw_dens$filenames

# Calculate the density of live cells
cleansummary = cleansummary[,dens_clean := (count_clean / count_unclean) * dens_unclean]

# Enforce date format
cleansummary$datetime_orig <- ymd_hms(cleansummary$datetime_orig)

# Write estimated cell densities
write.csv(cleansummary, './Script 6. Calculating cell density/output/density_dynamics.csv',
          row.names = FALSE)

# Remove all stored variables
rm(list = ls())
