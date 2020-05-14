rm(list = ls())
############################################################################################################
# SHANNON INFORMATION PARTITIONS: by Martin Hilbert & Mahima Agarwal @UCDavis.edu (unreleased version 0.0) #
############################################################################################################
#  It was adopted from an earlier Python implementation by Ryan James: https://github.com/Autoplectic/dit 
############################################################################################################
# Imports

# change the working directly to yours HERE
setwd("C:/Users/hilbert/OneDrive - University of California, Davis/Analytics Software/R_multivariate/ShannonPartition/code")
#File imports
source('distribution.R')
source('functions.R')
source('partitions.R')

#Package imports
# Note: You may need to install these packages first.
#   Do this using the install.packages('...') command 
library(igraph)
library(pander)
library(plyr)
library(purrr)
library(stringr)
###############################################################################
# Define constants


RV.MODES <- list(INDICES='index', NAMES='name')
history<- 1
###############################################################################
# Functions

# This function calculates the frequency of combinations from ts
counts.from.data<- function(ts, l)
{
  reshaped.ts<- matrix(NA, ncol=l, nrow=length(ts)-l+1)
  for (i in 1:length(ts)-l+1)
  {
    reshaped.ts[i,]<- c(ts[i:(i+l-1)])
  }
  reshaped.ts<- as.data.frame(reshaped.ts)
  counts<- count(reshaped.ts, vars = c('V1', 'V2'))
  return(counts)
}

###############################################################################
# OPTIONS 1: Input data manually
# Specify distribution of all possible joint realizations: you can input the data directly here. 
# Just specify the different combinations of variables that appear together, and then (in the same sequence) the joint probability of each (which should sum up to 1)
d<- create.distribution(matrix(c('WHA','WLA','WLB','SHA','SHC','SLD',
                                0.25,0.15,0.1,0.05,0.2,0.25), ncol=2, byrow=FALSE))


# Print all Shannon Information Partitions. THe different variables will be numbered 1,2,3... where 1 referes to the first variable you listed above, 2 the second, etc.
print.table(Shannon.Partition(d))

###############################################################################
# OPTIONS 2: Load two data time series
# have a folder named "data" inside the same folder where this script is, where the csv data file is
# change the name of your input file HERE 
data<- read.csv('../data/CatholicChurch.csv')
data.1<- data$lst_A.7
data.2<- data$lst_A.17

#This function creates the input distribution from the loaded csv data
distribution.from.data<- function(ts, l)
{
  counts<- counts.from.data(ts, l)
  counts.sum<- sum(counts$freq)
  pmf<- counts$freq/counts.sum
  create.distribution(outcomes=counts[,1:2], pmf=pmf) 
}

# Create a distribution from the data 
ts<- apply(cbind(data.1, data.2), 1, function(x){paste('(', paste(x, collapse=', '), ')', sep='')})
dist<- distribution.from.data(ts, history+1)

# Improve the distribution format by updating format of outcomes
dist<- distribution.modify.outcomes(dist, function(x){c(paste('(', get.tuple.nos(x[1]), ',)', sep=''),
                                                        get.tuple.nos(x[2]))})
dist<- distribution.set.rv.names(dist, c('X_past', 'Y_past', 'X_present', 'Y_present'))

# Print all Shannon information partitions
print.table(Shannon.Partition(dist), 7)

#Fetch specific subsets of partitions
get.entropy.subset(Shannon.Partition(dist), 'X_present', 'X_past')
get.entropy.subset(Shannon.Partition(dist), 'Y_past', 'X_past')
get.entropy.subset(Shannon.Partition(dist), 'Y_present', 'X_past')

###############################################################################
# OPTION 3: Load any number of time series (generalized case)
# change the name of your input file HERE
data<- read.csv('../data/CatholicChurch.csv')

# Create a distribution from the data 
ts<- apply(data, 1, function(x){paste('(', paste(x, collapse=', '), ')', sep='')})
dist<- distribution.from.data(ts, history + 1)

# Improve the distribution format by updating format of outcomes
dist<- distribution.modify.outcomes(dist, function(x){c(paste('(', get.tuple.nos(x[1]), ',)', sep=''),
                                                        get.tuple.nos(x[2]))})

colnames<- c(paste(colnames(data), 'past', sep='_'), paste(colnames(data), 'present', sep='_'))
dist<- distribution.set.rv.names(dist, colnames)

# Print all Shannon information partitions
partitions<- Shannon.Partition(dist) 
print.table(partitions)

#Fetch specific subsets of partitions: Examples
#Specify which aggregate of different infomration atoms you want to have printed out HERE 
get.entropy.subset(partitions, colnames[1], colnames[3])
get.entropy.subset(partitions, colnames[4], colnames[3])
get.entropy.subset(partitions, colnames[1], colnames[3])

###############################################################################