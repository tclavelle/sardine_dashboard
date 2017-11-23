##################################################
## Project: Sardine Dashboard V2
## Script purpose: Global script for shiny app
## Date: 11/22/2017
## Author: Tyler Clavelle
##################################################


# Packages ----------------------------------------------------------------

library(tidyverse)
library(dygraphs)
library(shiny)
library(scales)
library(TropFishR)
library(plotly)
library(xts)
library(lubridate)
library(knitr)
library(R.utils)
library(rsoi)

lwA <- 0.0078 # alpha from length-weight relationship
lwB <- 3.165 # beta from length weight relationship
t0 <- -lwA/lwB
max_age <- 31
recruit_age <- 7
ages <- c(7:max_age)
