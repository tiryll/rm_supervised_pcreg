# packages
                                  
install.packages(setdiff(c("docstring", "mvtnorm", "superpc", "pls",
                           "xtable","tikzDevice","ggplot2", "dplyr", "tidyr","readr"),
                         rownames(installed.packages())))

library(docstring)
library(mvtnorm)
library(superpc)
library(pls)
library(xtable)
library(tikzDevice)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# paths

data_dir <- paste(getwd(),"/bld/data",sep="")

tables_dir <- paste(getwd(),"/bld/tables",sep="")

figures_dir <- paste(getwd(),"/bld/figures",sep="")