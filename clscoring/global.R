library(shiny)
library(randomForest)
library(ROCR)
library(knitr)


source("rf_functions_toimport.R")

load("for_app.Rdata")


set.seed(2001)