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


### Script 5. Create a subset with data from all cleaned files equally represented
#     (this new subset will then be used for clustering to identify functional groups)
# Code takes as input:
#   1) Cleaned CSV files containing SFCM measurements (of live cells only), from Script 4
# See associated manuscript for details

#Set working directory
setwd()

#Set seed for reproducible results
set.seed(42)

#Load packages
library(data.table)
library(lubridate)

#Load list of filenames
filenames <- list.files(path = './Script 4. Cleaning by clustering/output/', 
                        pattern = c('cleaned'))

#Make empty data table to populate with values from each data file
numberoffiles <- length(filenames)
numberofpoints <- 100000L
pointsperfile <- as.integer(numberofpoints / numberoffiles)
numberofparameters <- 95L + 2L

m <- matrix(1, nrow = (numberoffiles * pointsperfile), ncol = numberofparameters)
dt <- as.data.table(m)
rm(m)

# Make the final row (datetime) a POSIXct object (origin present just for syntax) to avoid warning messages
dt <- dt[, V97 := as.POSIXct(V97, origin = '1950-06-27 13:00')]

# Create a list of row IDs where each new file's values will start at
idseq <- seq(1,(numberoffiles * pointsperfile), pointsperfile)

# Loop through all files to extract the data and populate the data table dt
for (index in 1:numberoffiles) {
  print(index)
  
  # Read data file. fread is a fast read function in the data.table package. 
  dat <- fread(paste0('./Script 4. Cleaning by clustering/output/', filenames[index]), 
               drop = c('id', 'biovol_Red2Total'))

  # Take a random subset from the file being read
  dat2 = dat[sample(.N, pointsperfile, replace = TRUE)]

  # Adding depth information (specific to this dataset)
  dat3 = dat2[, depth := as.numeric(as.character(substr(filenames[index], 32, 34)))]

  # Adding time information (specific to this dataset)
  dat4 = dat3[, datetime_orig := ymd_hm(gsub('u',':', substr(filenames[index], 15, 30)))]

  # populate existing data table with values
  set(dt, i = idseq[index]:(idseq[index] + pointsperfile - 1), 
      j = 1:(numberofparameters), dat4)
}

#Rename the (uninformative) column names in dt to those from the last data file
setnames(dt, names(dt), names(dat4))

#Write data subset
write.csv(dt, './Script 5. Generating cleaned data subset/output/clean_subset_100k.csv', 
          row.names = FALSE)

#Remove all objects from memory
rm(list = ls())
