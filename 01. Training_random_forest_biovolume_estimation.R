## Quantifying cell densities and biovolumes of phytoplankton communities and functional 
## groups using scanning flow cytometry, machine learning and unsupervised clustering
## Mridul K. Thomas, Simone Fontana, Marta Reyes & Francesco Pomati
#
## Code designed to be used with the dataset available at 10.5281/zenodo.977773
#
## Install all required packages before proceeding

list.of.packages <- c('data.table','randomForest', 'plot3Drgl', 'dplyr',
'lubridate', 'dtplyr', 'corrplot')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)

## This package has a different source and needs to be installed separately
source("http://bioconductor.org/biocLite.R")
biocLite("flowPeaks")


### Script 1. Training random forest to predict biovolume based on lab culture measurements
# Code takes as input:
#   1)  SFCM training dataset generated using lab cultures. Lab cultures were measured by 
#       SFCM and also by microscopy (for cell biovolumes). Non-cell particles were removed
#       from the SFCM data (i.e. the data was cleaned) by visually distinguishing cells 
#       from non-cell particles
# See associated manuscript for details

#set working directory
setwd()

# Load required packages
library(randomForest)

# Read in training dataset
culturedat <- read.csv('./Script 1. Training biovolume estimation/input/Random_forests_cleaned_lab_culture_trainingdata_2016_11_07.csv')

# Define new traits
culturedat$Red1Red2.ratio <- culturedat$FL.Red.Range / culturedat$X2.FL.Red.Range
culturedat$FL.Red.Gradient <- abs(culturedat$FL.Red.First - culturedat$FL.Red.Last) + 0.1
culturedat$X2.FL.Red.Gradient <- abs(culturedat$X2.FL.Red.First - culturedat$X2.FL.Red.Last) + 0.1

# Remove columns that do not capture information about cell physiology
culturedat$id <- NULL
culturedat$SWS.Time.of.Arrival <- NULL
culturedat$SWS.Sample.Length <- NULL

# Train random forest to estimate biovolume based on all available traits
rf <- randomForest(culturedat[,5:98], culturedat$microscopy.biovolume,
                       importance = TRUE, ntree = 1001)

saveRDS(rf, './Script 1. Training biovolume estimation/output/Trained_RF_Cytobuoy_biovolume_estimation.rds')
varImpPlot(rf, type = 1)

# Remove all stored variables
rm(list = ls())
