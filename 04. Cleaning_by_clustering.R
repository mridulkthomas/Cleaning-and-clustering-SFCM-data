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


### Script 4. Cluster raw data subset using variables identified in Script 2. 
#     (Clusters data, then writes clean data files without points belonging to clusters 
#     visually identified as 'junk' i.e. not living phytoplankton cells. Also estimates
#     biovolumes of individual cells)
# Code takes as input: 
#   1) Raw data files containing SFCM measurements
#   2) Raw data subset, from Script 3
#   3) Variables identified as being important in distinguishing live cells from 
#       non-cell particles, from Script 2
#   4) Random forest for biovolume estimation, from Script 1
# See associated manuscript for details

#Load packages

library(flowPeaks)
library(data.table)
library(randomForest)
library(plot3Drgl)

#Set seed for repeatability
set.seed(42)

#Set your working directory, which contains both the data subset and the files to be cleaned (the latter can be in subfolders)
setwd()

#Define colour palette for clusters in 3D plots
palette1 <- c('black','red','green3','blue','cyan','magenta','yellow','gray',
              'brown','coral4','khaki','orange','darkgreen')
palette(palette1)

# Using as input the output from Scipt 1 and Script 3.
# Read in stored random forests
rf_biovol <- readRDS('./Script 1. Training biovolume estimation/output/Trained_RF_Cytobuoy_biovolume_estimation.rds')

# Read in the data subset to be used for clustering. Change filename as needed.
df <- read.csv('./Script 3. Generating raw data subset/output/raw_data_subset_100k.csv')

# Remove columns that are believed to be inaccurate
df$Length.Max <- df$Length.Min <- df$Length.Max.Corrected <- NULL

# Define new variables to help with identification of electonic noise
# 0.1 added to all Gradient values because a) there are zeroes, b) it is just below the lowest non-zero measurement, c) we need to log transform later
df$FL.Red.Gradient <- abs(df$FL.Red.First - df$FL.Red.Last) + 0.1
df$X2.FL.Red.Gradient <- abs(df$X2.FL.Red.First - df$X2.FL.Red.Last) + 0.1
df$Red1Red2.ratio <- df$FL.Red.Range / df$X2.FL.Red.Range

# Select only the variables used for clustering and visualization
dat.clust <- df[,c('FWS.Range','FL.Red.Range','X2.FL.Red.Range', 
                   'FL.Yellow.Range', 'X2.FL.Red.Gradient',
                   'FL.Orange.Range', 'FL.Red.Number.of.cells', 
                   'X2.FL.Red.Last', 'FL.Red.Fill.factor',
                   'X2.FL.Red.Fill.factor')]

# Subset data to remove extreme cases that are clearly not live cells 
# i.e. below size range of phytoplankton cells or showing near-zero pigment fluorescence.
# Extreme outliers will influence  clustering.
dat.clust <- subset(dat.clust,dat.clust$FL.Red.Range > 0.2)
dat.clust <- subset(dat.clust,dat.clust$X2.FL.Red.Range > 0)
dat.clust <- subset(dat.clust,dat.clust$FL.Yellow.Range > 0)
dat.clust <- subset(dat.clust,dat.clust$FL.Orange.Range > 0)

# Exclude additional cases based on plots and individual judgement.
# Here we exclude a few extreme outliers which negatively influence clustering results. 
# Visual examination of plots indicates that these are not cells
dat.clust <- subset(dat.clust, dat.clust$FWS.Range > 0.2)

# Log-transform all variables
dat.clustlog <- log10(dat.clust)
# =================== MSC ==================
# ============== TODO @Nik @Kel ============
# realize how they determine bad-good clusters (e.g high fluorescence on one channel and very low on others)
# instead of performing clustering here, run the whole dataset through the random forest calculated in step 2 
#           (so we have to load the whole dataset and not just the raw dataset subset from script 3
# if that takes a lot of time, re-construct the random forest using the "most-important parameters", those close to the root node
# Compare with auto-encoder
# save this "cleaned dataset"
# after removing all those points, then perform clustering on a subset from the "cleaned" raw dataset 
#            - so we have to run Script 3 after this step 
# visualize the clusters and validate with training set (classify training set on clusters, measure accuracy)
# investigate other clustering methods other than flowPeaks
# ==========================================
# ==========================================
# Cluster using flowPeaks

subflow <- flowPeaks(dat.clustlog[,c(1:10)], tol = 0.25, h0 = 0.05, h = 2) 


# Examine details if interested. 
# 'weight' column indicates the proportion of the data belonging to that cluster
summary(subflow)

# Examine 3D plots to 1) verify whether clusters are meaningful and 2) identify clusters
# corressponding to live phytoplankton cells
# Junk clusters  are characterized by low FWS.Range and low fluorescence.
# High fluorescence on one channel but very low fluorescence on other channels can indicate electronic noise.
# Typically the largest cluster is junk, although there may be more than one

scatter3Drgl(dat.clustlog$FWS.Range, dat.clustlog$FL.Red.Range, dat.clustlog$X2.FL.Red.Range,
             colvar = subflow$peaks.cluster, col = 1:max(subflow$peaks.cluster),
             xlab = 'FWS.Range', ylab = 'FL.Red.Range', zlab = 'X2.FL.Red.Range')

scatter3Drgl(dat.clustlog$FL.Red.Range, dat.clustlog$X2.FL.Red.Range, dat.clustlog$X2.FL.Red.Gradient,
             colvar = subflow$peaks.cluster, col = 1:max(subflow$peaks.cluster),
             xlab = 'FL.Red.Range', ylab = 'X2.FL.Red.Range', zlab = 'X2.FL.Red.Gradient')

# Enter the cluster numbers of the junk clusters below (based on 3D plots)
# For this dataset, cluster #1 alone appears to be junk
junk.clusterids <- c(1)
dat.clustlog2 <- data.frame(dat.clustlog,subflow$peaks.cluster)

# Create dataset excluding junk clusters
dat.good <- subset(dat.clustlog2, !(dat.clustlog2$subflow.peaks.cluster %in% junk.clusterids))

# Re-examine plot without the putative junk clusters (this may reveal smaller ones that were not visible earlier).
# Note that axis ranges & colours will not remain consistent with previous set of plots
scatter3Drgl(dat.good$FWS.Range, dat.good$FL.Red.Range, dat.good$X2.FL.Red.Range,
             colvar = dat.good$subflow.peaks.cluster, col = 1:max(subflow$peaks.cluster),
             xlab = 'FWS.Range', ylab = 'FL.Red.Range', zlab = 'X2.FL.Red.Range')

scatter3Drgl(dat.good$X2.FL.Red.Range, dat.good$FL.Red.Range, dat.good$X2.FL.Red.Gradient,
             colvar = dat.good$subflow.peaks.cluster, col = 1:max(subflow$peaks.cluster),
             xlab = 'X2.FL.Red.Range', ylab = 'FL.Red.Range', zlab = 'X2.FL.Red.Gradient')

# When satisfied with clustering results, train a random forest classifier to distinguish the clusters
rf_cleaning <- randomForest(factor(subflow.peaks.cluster) ~ FWS.Range + FL.Red.Range + 
                              X2.FL.Red.Range + FL.Yellow.Range + FL.Orange.Range + 
                              X2.FL.Red.Gradient + FL.Red.Number.of.cells + X2.FL.Red.Last +
                              FL.Red.Fill.factor + X2.FL.Red.Fill.factor, 
                            data = dat.clustlog2,importance = TRUE, ntree = 1001, 
                            na.action = na.omit)

# Remove objects that are no longer needed
rm(subflow, dat.clust, dat.clustlog, dat.clustlog2, dat.good, df, palette1)

# Load the list of raw data files
filenames <- list.files(path = './Script 3. Generating raw data subset/input/', 
                        pattern = c('Allparameters_'), full.names = FALSE)

# Define function to write cleaned files (i.e. excluding junk particles) 
# using trained random forest, and to estimate biovolumes of cells in cleaned files.

classifydat <- function(x){
  # Read in data file, exclude unhelpful columns unless needed
  dat <- fread(paste0('./Script 3. Generating raw data subset/input/', x), 
               drop = c('V1','SWS.Sample.Length', 'SWS.Time.of.Arrival', 'Length.Max',
                           'Length.Min','Length.Max.Corrected'))

  # Impose same conditions as at the beginning of this process, or clustering may go awry.
  dat = subset(dat,FL.Red.Range > 0.2)
  dat = subset(dat,X2.FL.Red.Range > 0)
  dat = subset(dat,FL.Orange.Range > 0)
  dat = subset(dat,FL.Yellow.Range > 0)
  dat = subset(dat,FWS.Range > 0.2)
  
  # Define new parameters
  dat = dat[, FL.Red.Gradient := abs(FL.Red.First - FL.Red.Last) + 0.1]
  dat = dat[, X2.FL.Red.Gradient := abs(X2.FL.Red.First - X2.FL.Red.Last) + 0.1]
  dat = dat[, Red1Red2.ratio := (FL.Red.Range / X2.FL.Red.Range)]

  # Arrange by id.
  setkey(dat, id)

  # Select only the variables used for clustering and log-transform
  dat2 <- dat[,list(FWS.Range, FL.Red.Range, X2.FL.Red.Range, FL.Yellow.Range, 
                   X2.FL.Red.Gradient, FL.Orange.Range, FL.Red.Number.of.cells, 
                   X2.FL.Red.Last, FL.Red.Fill.factor, X2.FL.Red.Fill.factor)]
  dat2 = log10(dat2)

  # Predict the cluster identity of particles using the random forests classifier
  preds_clean <- predict(rf_cleaning, dat2)

  # Add a new column with the predicted cluster identity
  dat[,classification1 := as.numeric(preds_clean)]

  # Remove the clusters that are believed to be junk
  dat3 = dat[!(classification1 %in% junk.clusterids)]

  # Remove classification columns if desired.
  dat3[, classification1 := NULL]

  # Predict biovolumes of every cell with a random forest
  preds_biovol <- predict(rf_biovol, dat3)

  # Add biovolume predictions to dataset
  datfinal = dat3[, biovol_rf := preds_biovol]

  # Arrange by id (i.e. original measurement sequence). 
  # Not necessary, but helpful for later error checking
  setkey(datfinal, id)

  # Write cleaned file to .csv
  # Cleaned files are written to same subfolder as original file
  write.csv(datfinal, paste0('./Script 4. Cleaning by clustering/output/', 
                             substr(x, 1, nchar(x) - 4), '_cleaned', '.csv'),
            row.names = FALSE)

}

# Run function on file list to clean sequentially.
mapply(classifydat, x = filenames)

# Remove all stored variables
rm(list = ls())
