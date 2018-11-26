library(shiny)

library(randomForest)
#library(reprtree)
library(ROCR)
library(knitr)
#library(maftools)

#source("C:/Users/Vladimir/Desktop/tcga/rf_functions_toimport.R")

source("rf_functions_toimport.R")

#load("../data/for_app.Rdata")

load("for_app.Rdata")


set.seed(2001)