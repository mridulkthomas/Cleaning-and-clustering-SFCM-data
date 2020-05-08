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

###kel-In this dataset there are the complete records of the field dataset, the samples were taken at 6 depths every 4 hours
#for 71 days. So we have totally 6x71 = 426 SFCM measurements##


### Script 3. Create a subset with data from all raw files equally represented
#     (this subset will then be used for clustering and then to remove junk particles)
# Code takes as input: 
#   1) Raw data files containing SFCM measurements
# See associated manuscript for details

#set working directory
setwd()

#Set seed for reproducible clustering results
set.seed(42)

#Load required packages
library(data.table)

#Load list of filenames for subsetting
filenames <- list.files(path = 'Script 3. Generating raw data subset/input/', 
                        pattern = c('Allparameters_'), recursive = TRUE)

#Make empty data table to populate with values from each data file
numberoffiles <- length(filenames)
numberofpoints <- 100000L
pointsperfile <- as.integer(numberofpoints / numberoffiles)
numberofparameters <- 96L

m <- matrix(1,nrow = (numberoffiles * pointsperfile), ncol = numberofparameters)
dt <- as.data.table(m)

#Make the first row (id) an integer to avoid warning messages
dt <- dt[, V1 := as.integer(V1)]

#Removing the matrix
rm(m)

#Create a list of row IDs where each new file's values will start at
idseq <- seq(1,(numberoffiles * pointsperfile), pointsperfile)

#Loop through all files to extract the data and populate the data table dt

for (index in 1:numberoffiles) {
  print(index)
  # Read data file. fread is a fast read function in the data.table package. 
  dat <- fread(paste0('Script 3. Generating raw data subset/input/', 
                      filenames[index]), drop = c('V1', 'id'))
  
  # Remove rare cases where fluorescence range = 0, since they are clearly 
  # not live cells and we will log transform the data later
  dat2 = dat[dat$FL.Red.Range > 0]

  # Take a random subset from the file being read
  dat3 = dat2[sample(.N, pointsperfile)]

  # Populate existing data table with values
  set(dt, i = idseq[index]:(idseq[index] + pointsperfile - 1), 
      j = 1:(numberofparameters), dat3)
}

rm(index)

#Renames the (uninformative) column names in dt to those from the last data file
setnames(dt, names(dt), names(dat))

# Write data subset
write.csv(dt, './Script 3. Generating raw data subset/output/raw_data_subset_100k.csv', 
          row.names = FALSE)

# Remove all stored variables
rm(list = ls())
