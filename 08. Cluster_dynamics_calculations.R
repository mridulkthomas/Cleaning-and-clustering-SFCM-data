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


### Script 8. Cluster cleaned data subset and estimate the cell densities and total 
#       biovolumes of individual clusters in the cleaned data files.
# Code takes as input: 
#   1) Cleaned data files, from Script 4
#   2) Cleaned data subset, from Script 5
#   3) Variables identified as being important in distinguishing functional groups, from
#       Script 7
# See associated manuscript for details

#Set working directory
setwd()

#Set seed for repeatability
set.seed(42)

#Load packages
library(flowPeaks)
library(lubridate)
library(data.table)
library(randomForest)
library(plot3Drgl)
library(dplyr)
library(dtplyr)

# Create palette for cluster plots
palette1 <- c('black','red','green3','blue','cyan','magenta','yellow','gray','brown',
              'purple', 'goldenrod','pink','orange','darkgreen', 'coral2', 'khaki4')
palette(palette1)

# Read in the data subset to be used for clustering
df <- read.csv('./Script 5. Generating cleaned data subset/output/clean_subset_100k.csv')

# Exclude cells too long or short for the instrument to measure; indicative of a system error
df <- df[df$FWS.Length < 1000 & df$FWS.Length > 2,]

# Select only the variables used for clustering
dat.clust <- df[, c('FL.Red.Range','X2.FL.Red.Range','FL.Orange.Range',
                    'Red1Red2.ratio', 'X2.FL.Red.Gradient',
                    'FWS.Length','FWS.Fill.factor', 'FL.Red.First')]

# Take the log of all except for FWS.Fill.factor, which ranges from 0 to 1
dat.clustlog <- log10(dat.clust)

# Cluster using flowPeaks
subflow <- flowPeaks(dat.clustlog[, c(1:8)], tol = 0.25, h0 = 0.05, h = 2)

# Examine details if interested
summary(subflow)

# Examine any plots of interest
scatter3Drgl(dat.clustlog$Red1Red2.ratio, dat.clustlog$FWS.Length, 
             dat.clustlog$FL.Orange.Range,
             colvar = subflow$peaks.cluster, col = 1:length(unique(subflow$peaks.cluster)),
             xlab = 'Red1/Red2 ratio', ylab = 'FWS.Length', zlab = 'FL.Orange.Range')

scatter3Drgl(dat.clustlog$Red1Red2.ratio, dat.clustlog$X2.FL.Red.Gradient, 
             dat.clustlog$FWS.Fill.factor,
             colvar = subflow$peaks.cluster, col = 1:length(unique(subflow$peaks.cluster)),
             xlab = 'Red1/Red2 ratio', ylab = 'X2.FL.Red.Gradient', zlab = 'FWS.Fill.factor')

# Add cluster identity to dataset
dat.clustlog$subflow.peaks.cluster <- subflow$peaks.cluster

# Train a random forest classifier to identify the clusters
rf_groups <- randomForest(factor(subflow.peaks.cluster) ~ FWS.Length + FL.Red.Range + 
                            X2.FL.Red.Range + FL.Orange.Range + Red1Red2.ratio + 
                            X2.FL.Red.Gradient + FWS.Fill.factor + FL.Red.First, 
                          data = dat.clustlog,importance = TRUE, ntree = 1001, 
                          na.action = na.omit)

# Load the list of cleaned data files
cleanfiles <- list.files(path = './Script 4. Cleaning by clustering/output/', 
                         pattern = c('cleaned'), recursive = T, full.names = TRUE)

# Generate empty file to be populated with summary data
numberoffiles <- length(cleanfiles)
m <- matrix(1, nrow = (numberoffiles * length(rf_groups$classes)), ncol = 5)
dt <- as.data.table(m)
rm(m)

# Number of clusters identified
cluster_number <- length(rf_groups$classes)

# Fill in cluster number column
set(dt, i = NULL, j = 1L, rep(1L:cluster_number,length(cleanfiles)))

# Fill in column names
setnames(dt, names(dt), c('cluster','cluster_count','cluster_percent', 'biovol_rf_mean',
                          'datetime_orig'))

# Set non-numeric column types
dt <- dt[, cluster_count := as.integer(cluster_count)]
dt <- dt[, datetime_orig := as.character(datetime_orig)]

# Create list of starting row values for each file (i.e. start after every cluster has a row)
idseq <- seq(1,(numberoffiles * cluster_number), cluster_number)

# Create a list of all clusters for merging later (otherwise missing clusters are omitted)
clust_list <- data.frame(cluster = seq(1:cluster_number))

#Loop through all files
for (index in 1:numberoffiles) {

  print(index)

  #Read in data file. Only read columns for group clustering and for biovolume calculations
  dat <- fread(cleanfiles[index], select = c('FL.Red.Range','X2.FL.Red.Range',
                                             'FL.Orange.Range', 'Red1Red2.ratio', 
                                             'X2.FL.Red.Gradient', 'FWS.Length',
                                             'FWS.Fill.factor', 'FL.Red.First',
                                             'biovol_rf'))

  # Log-transform relevant variables and add a new one to classify with random forest
  dat2 <- log10(dat)

  # Predict the cluster identity of particles using the random forests classifier
  preds <- predict(rf_groups, dat2)

  # Add cluster IDs to dataset
  dat3 = dat[, cluster := preds]

  # Calculate biovolume by group
  biovol_rf_mean <- group_by(dat3, cluster) %>%
                      filter(FWS.Length > 2, FWS.Length < 1000) %>%
                      summarise(mean(biovol_rf)) %>%
                      merge(., clust_list, all.y = TRUE, by = 'cluster') %>%
                      rename(value = `mean(biovol_rf)`)

  # Save number of cells belonging to each cluster
  set(dt, i = idseq[index]:(idseq[index] + cluster_number - 1), j = 2L, summary(preds))

  # Save proportion of cells belonging to each cluster
  set(dt, i = idseq[index]:(idseq[index] + cluster_number - 1), j = 3L,
      summary(preds) / sum(summary(preds)))

  # Save mean biovolume (by random forest) of each cluster
  set(dt, i = idseq[index]:(idseq[index] + cluster_number - 1), j = 4L, biovol_rf_mean$value)

  #Save datetime information
  set(dt, i = idseq[index]:(idseq[index] + cluster_number - 1), j = 5L,
      rep(as.character(ymd_hm(gsub('u',':', substr(cleanfiles[index], 57, 72)))),
          cluster_number))

}

# Remove unnecessary variables and data frames
rm(index, clust_list, cleanfiles, cluster_number, idseq, numberoffiles, rf_groups, 
   dat, dat2, dat3, biovol_rf_mean)

# Replace NAs (missing values) with zeroes for total biovolume columns
dt$biovol_rf_mean[is.na(dt$biovol_rf_mean)] <- 0

# Enforce date format 
dt$datetime_orig <- ymd_hms(dt$datetime_orig)

# Read in file with density calculations
densdat <- read.csv('./Script 6. Calculating cell density/output/density_dynamics.csv')
densdat$datetime_orig <- dmy_hm(densdat$datetime_orig) #might need to change to ymd_hms

# Merge cluster and density files for calculation of density and biovolume percentages of each cluster
clustdat <- merge(dt, densdat, by = 'datetime_orig')
clustdat$dens_clust <- clustdat$dens_clean * clustdat$cluster_percent
clustdat$totalbiovol_rf <- clustdat$biovol_rf_mean * clustdat$dens_clust

# Write file
write.csv(clustdat, './Script 8. Calculating cell densities of clusters/output/cluster_dynamics.csv', row.names = FALSE)

# Remove all stored variables
rm(list = ls())
