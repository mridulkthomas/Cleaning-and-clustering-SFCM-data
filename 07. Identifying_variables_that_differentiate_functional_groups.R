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


### Script 7. Identify variables that help to distinguish functional groups from each other
#     (The top variables will then be used for clustering to identify clusters within the
#       cleaned SFCM data)
# Code takes as input:
#   1) SFCM training dataset generated using lab cultures (previously used in Script 1 
#         for biovolume estimation). Lab cultures were measured by SFCM and non-cell
#         particles were removed from the SFCM data (i.e. the data was cleaned) by visually
#         distinguishing cells from non-cell particles
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
dat_init <- read.csv('./Script 1. Training biovolume estimation/input/Random_forests_cleaned_lab_culture_trainingdata_2016_11_07.csv')

# Defining new variables
dat_init$FL.Red.Gradient <- abs(dat_init$FL.Red.First - dat_init$FL.Red.Last) + 0.1
dat_init$X2.FL.Red.Gradient <- abs(dat_init$X2.FL.Red.First - dat_init$X2.FL.Red.Last) + 0.1
dat_init$Red1Red2.ratio <- dat_init$FL.Red.Range / dat_init$X2.FL.Red.Range

# Remove columns that do not contain SFCM information about cell physiology
dat_init$id <- NULL
dat_init$SWS.Time.of.Arrival <- NULL
dat_init$SWS.Sample.Length <- NULL
dat_init$microscopy.biovolume <- NULL
dat_init$microscopy.length <- NULL

# Log-transforming variables without zeroes or negative values
dat_logged <- log10(dat_init[, -grep("Curvature|Asymmetry|species|group", names(dat_init))])
dat_unlogged <- dat_init[, grep("Curvature|Asymmetry|species|group", names(dat_init))]

dat_allvars <- cbind(dat_unlogged, dat_logged)
rm(dat_logged, dat_unlogged)

# Identifying sets of parameters that are highly correlated 
cormat <- cor(dat_allvars[,3:96])
corrplot(cormat,  type = 'lower', method = 'number', addCoefasPercent = TRUE, 
         tl.cex = 0.5, number.cex = 0.5, diag = FALSE)

# Removing highly correlated variables
dat <- dat_allvars[, -grep("Total|Average|Maximum|Inertia|Center.of.gravity", 
                                names(dat_allvars))]

dat$SWS.Length <- NULL
dat$SWS.Number.of.cells <- NULL
dat$FL.Red.Number.of.cells <- NULL
dat$FL.Red.Length <- NULL
dat$FL.Yellow.Length <- NULL
dat$SWS.Fill.factor <- NULL
dat$FWS.First <- NULL
dat$SWS.First <- NULL
dat$FL.Orange.First <- NULL
dat$FL.Red.Last <- NULL
dat$FL.Yellow.First <- NULL
dat$SWS.Range <- NULL
dat$X2.FL.Red.Last <- NULL
dat$X2.FL.Red.Asymmetry <- NULL
dat$FL.Yellow.Range <- NULL
dat$FL.Orange.Last <- NULL

# Training random forest to identify most important variables in distinguishing between
# phytoplankton functional groups
rf <- randomForest(dat[,3:45], dat$group, importance = TRUE, ntree = 10001)
varImpPlot(rf, type = 2)

saveRDS(rf, 'E:/Cytobuoy_data/RF_lab_measurements_distinguishing_groups_variable_importance.rds')

# Saving trained random forest
saveRDS(rf, './Script 7. Identifying important variables for group distinction/output/Trained_RF_lab_measurements_data_cleaning_variable_importance.rds')

# Remove all stored variables
rm(list = ls())
