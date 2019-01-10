library(shiny)
library(randomForest)
library(ROCR)
library(knitr)
library(tidyverse); library(magrittr); library(purrr)
library(png); library(grid)

#load("for_app.Rdata")
load("for_app_BRCA.Rdata")
load("for_app_UCEC.Rdata")
proc_img <- readPNG("img/processing.png")

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