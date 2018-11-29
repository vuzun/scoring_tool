library(shiny)
library(randomForest)
library(ROCR)
library(knitr)
library(tidyverse); library(magrittr); library(purrr)

#load("for_app.Rdata")
load("for_app_BRCA.Rdata")
load("for_app_UCEC.Rdata")

source("rf_functions_toimport.R")


set.seed(2001)