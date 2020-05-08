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

#library(flowPeaks)
library(data.table)
library(randomForest)

#Set seed for repeatability
set.seed(42)
#set working directory
setwd("/cloud/project")

#Define colour palette for clusters in 3D plots
palette1 <- c('black','red','green3','blue','cyan','magenta','yellow','gray',
              'brown','coral4','khaki','orange','darkgreen')
palette(palette1)


# Load the list of raw data files
filenames <- list.files(path = './Script 3. Generating raw data subset/input/', 
                        pattern = c('Allparameters_'), full.names = FALSE)

clean_dat <- function(x){
  # Read in data file, exclude unhelpful columns unless needed
  df <- fread(paste0('./Script 3. Generating raw data subset/input/', x), 
              drop = c('V1','SWS.Sample.Length', 'SWS.Time.of.Arrival', 'Length.Max',
                       'Length.Min','Length.Max.Corrected'))
  
  # Remove columns that are believed to be inaccurate
  df$Length.Max <- df$Length.Min <- df$Length.Max.Corrected <- NULL
  
  # Define new variables to help with identification of electonic noise
  # 0.1 added to all Gradient values because a) there are zeroes, b) it is just below the lowest non-zero measurement, c) we need to log transform later
  df$FL.Red.Gradient <- abs(df$FL.Red.First - df$FL.Red.Last) + 0.1
  df$X2.FL.Red.Gradient <- abs(df$X2.FL.Red.First - df$X2.FL.Red.Last) + 0.1
  df$Red1Red2.ratio <- df$FL.Red.Range / df$X2.FL.Red.Range
  
  # Select only the variables used for clustering and visualization
  df <- df[,c('FWS.Range','FL.Red.Range','X2.FL.Red.Range', 
              'FL.Yellow.Range', 'X2.FL.Red.Gradient',
              'FL.Orange.Range', 'FL.Red.Number.of.cells', 
              'X2.FL.Red.Last', 'FL.Red.Fill.factor',
              'X2.FL.Red.Fill.factor')]
  
  # Subset data to remove extreme cases that are clearly not live cells 
  # i.e. below size range of phytoplankton cells or showing near-zero pigment fluorescence.
  # Extreme outliers will influence  clustering.
  df <- subset(df,df$FL.Red.Range > 0.2)
  df <- subset(df,df$X2.FL.Red.Range > 0)
  df <- subset(df,df$FL.Yellow.Range > 0)
  df <- subset(df,df$FL.Orange.Range > 0)
  
  # Exclude additional cases based on plots and individual judgement.
  # Here we exclude a few extreme outliers which negatively influence clustering results. 
  # Visual examination of plots indicates that these are not cells
  df <- subset(df, df$FWS.Range > 0.2)
  # Log-transform all variables
  dflog <- log10(df)
  
  
  # Write cleaned file to .csv
  # Cleaned files are written to same subfolder as original file
  write.csv(df, paste0('./Script 4. Cleaning by clustering/output/', 
                       substr(x, 1, nchar(x) - 4), '_alt_cleaned', '.csv'),
            row.names = FALSE)
  
}

# Run function on file list to clean sequentially.
mapply(clean_dat, x = filenames)
