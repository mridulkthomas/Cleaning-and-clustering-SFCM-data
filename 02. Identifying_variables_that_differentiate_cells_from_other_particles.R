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

#kel-This script takes as input both live and non-live cells and puts labels on every sample.#

### Script 2. Identify variables that help to distinguish live cells from other particles
#     (the top variables will then be used for clustering to remove junk particles)
# Code takes as input:
#   1) SFCM training dataset generated using lab cultures. Cells and non-cell particles 
#       were distinguished visually. Identity of individual particles (cell / non-cell)
#       are indicated in column 'qual' (quality).
#       'good' = live cell, 'bad' = non-cell particle
# See associated manuscript for details

#Set seed for repeatability
set.seed(42)

#set working directory
setwd()

# Load required packages
library(data.table)
library(randomForest)
library(corrplot)

# Read in training dataset
dat_init <- read.csv('./Script 2. Identifying important variables for data cleaning/input/Random_forests_data_cleaning_trainingdata_2016_10_30.csv')

# Defining new variables
dat_init$FL.Red.Gradient <- abs(dat_init$FL.Red.First - dat_init$FL.Red.Last) + 0.1
dat_init$X2.FL.Red.Gradient <- abs(dat_init$X2.FL.Red.First - dat_init$X2.FL.Red.Last) + 0.1

# Log-transforming variables without zeroes or negative values
dat_logged <- log10(dat_init[, -grep("Curvature|Asymmetry|qual|species|group", names(dat_init))])
dat_unlogged <- dat_init[, grep("Curvature|Asymmetry|qual|species|group", names(dat_init))]

dat <- cbind(dat_logged, dat_unlogged)
rm(dat_logged, dat_unlogged)

# Identifying sets of parameters that are highly correlated 
cormat <- cor(dat[,1:94])
corrplot(cormat,  type = 'lower', method = 'number', addCoefasPercent = TRUE, 
         tl.cex = 0.5, number.cex = 0.5, diag = FALSE)

# Removing highly correlated variables
dat <- dat[, -grep("Total|Average|Maximum|Inertia", names(dat))]
dat <- dat[, -grep("Center.of.gravity", names(dat))]

dat$SWS.Fill.factor <- NULL
dat$SWS.First <- NULL
dat$FWS.First <- NULL
dat$FWS.Last <- NULL
dat$FL.Orange.Last <- NULL
dat$FL.Orange.Asymmetry <- NULL
dat$Curvature.Minimum <- NULL
dat$X2.FL.Red.Asymmetry <- NULL
dat$SWS.Length <- NULL
dat$SWS.Range <- NULL

# Training random forest to identify most important variables in distinguishing between
# cells and other particles
rf <- randomForest(dat[, 1:49], dat$qual, importance = TRUE, ntree = 10001)
varImpPlot(rf, type = 2)

# Saving trained random forest
saveRDS(rf, './Script 2. Identifying important variables for data cleaning/output/Trained_RF_lab_measurements_data_cleaning_variable_importance.rds')

# Remove all stored variables
rm(list = ls())
