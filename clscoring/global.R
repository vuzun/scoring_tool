library(shiny)
library(randomForest)
library(ROCR)
library(knitr)
library(tidyverse); library(magrittr); library(purrr)
library(png); library(grid)

#load("for_app.Rdata")
load("for_app_BRCA.Rdata")
load("for_app_UCEC.Rdata")
proc_img <- readPNG("../data/processing.png")

source("rf_functions_toimport.R")


set.seed(2001)

#---todo---
#
#AUC
#combo
#mds cell line
#WORKING CUSTOM SCHEMES
#mds doesnt show

#grid.raster(proc_img)
#vals$rfmds <<- recordPlot()}

#having separate hists for each tumour type, so the update overrides only one of them
#+button to set the default back