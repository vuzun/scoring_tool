library(shiny)
library(randomForest)
library(ROCR)
library(knitr)
library(tidyverse)
library(magrittr)
library(png); library(grid)

load("for_app_BRCA.Rdata")
load("for_app_UCEC.Rdata")

source("rf_functions_toimport.R")


set.seed(2001)



#grid.raster(proc_img)
#vals$rfmds <<- recordPlot()}
#having separate hists for each tumour type, so the update overrides only one of them
#+button to set the default back