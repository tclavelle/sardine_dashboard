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


# Data --------------------------------------------------------------------

comm_df <- read_csv(file = 'data/C80PNVPPcm.csv', na = '...') # quarterly commercial landings data
muni_df <- read_csv(file = 'data/D10PNVMPcm.csv', na = '...') # quarterly municipal landings data

lwA <- 0.0078 # alpha from length-weight relationship
lwB <- 3.165 # beta from length weight relationship
t0 <- -lwA/lwB
max_age <- 31
recruit_age <- 7
ages <- c(7:max_age)

# Catch data processing for dygraph ---------------------------------------

catch_df <- comm_df %>%
  mutate(Fleet = 'Commercial') %>%
  bind_rows(muni_df %>%
              mutate(Fleet = 'Municipal')) %>%
  mutate(Day   = '01',
         Month = NA)

catch_df$Month[catch_df$Date == 'Quarter 1'] <- '04'
catch_df$Month[catch_df$Date == 'Quarter 2'] <- '07'
catch_df$Month[catch_df$Date == 'Quarter 3'] <- '10'
catch_df$Month[catch_df$Date == 'Quarter 4'] <- '12'
catch_df$Day[catch_df$Date == 'Quarter 4'] <- '31'


catch_df <- catch_df %>%
  mutate(Date  = paste(Year, Month, Day, sep = '-'),
         Catch = as.numeric(Catch),
         Date  = as_date(Date)) %>%
  select(Date, Catch, Fleet) %>%
  ungroup()

### Catch history plot code
plot1 <- catch_df %>%
  filter(is.na(Catch)==F) %>%
  spread(Fleet, Catch) %>%
  mutate(All  = Commercial + Municipal,
         Season = ifelse(Date > '2011-12-31' & Date < '2017-04-01' & grepl('-12', Date), 0, 1)) # indicate closed season

p1_comm <- xts(plot1$Commercial, order.by = plot1$Date)
p1_muni <- xts(plot1$Municipal, order.by = plot1$Date)
p1_all <- xts(plot1$All, order.by = plot1$Date)
season <- xts(plot1$Season, order.by = plot1$Date)
p1 <- cbind(p1_muni, p1_all)
