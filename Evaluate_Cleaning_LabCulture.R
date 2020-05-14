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
library(googledrive)
#Set seed for repeatability
set.seed(42)
#set working directory
setwd("/cloud/project")

#@Nik Retrieve 8GB Field Dataset from Google Drive 
#(Accesss Rights & Google Authorization Required)

#googledrive::drive_auth(use_oob = TRUE)
file = googledrive::drive_get("Phytoplankton_Data/workable.zip")
print(file)

googledrive::drive_download(file)

#Define colour palette for clusters in 3D plots
palette1 <- c('black','red','green3','blue','cyan','magenta','yellow','gray',
              'brown','coral4','khaki','orange','darkgreen')
palette(palette1)


# Load the list of raw data files
filenames <- list.files(path = './workable/', 
                        pattern = c('Random_forests'), full.names = FALSE)

clean_dat <- function(x){
  # Read in data file, exclude unhelpful columns unless needed
  df <- fread(paste0('./workable/', x))
  print(colnames(df))
  # Remove columns that are believed to be inaccurate
  #df$Length.Max <- df$Length.Min <- df$Length.Max.Corrected <- NULL
  
  # Define new variables to help with identification of electonic noise
  # 0.1 added to all Gradient values because a) there are zeroes, b) it is just below the lowest non-zero measurement, c) we need to log transform later
  df$FL.Red.Gradient <- abs(df$FL.Red.First - df$FL.Red.Last) + 0.1
  df$X2.FL.Red.Gradient <- abs(df$X2.FL.Red.First - df$X2.FL.Red.Last) + 0.1
  df$Red1Red2.ratio <- df$FL.Red.Range / df$X2.FL.Red.Range

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
  write.csv(df, paste0('./Cleaned_Field/',x),
            row.names = FALSE)
}

dir.create('./Cleaned_Field/',showWarnings = FALSE)
# Run function on file list to clean sequentially.
out = mapply(clean_dat, x = filenames)

for(x in filenames){
  df <- fread(paste0('./Cleaned_Field/',x ))
  total_cells = length(rownames(df))
  print(df$qual)
  qual = table(df$qual)
  dead_cells = qual[names(qual) == "bad"]
  FP = (dead_cells/total_cells)*100
  print(paste0('False Positive = ',FP))
  
}













