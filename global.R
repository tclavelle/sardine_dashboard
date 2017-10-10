##################################################
## Project: EDF Sardine Tool
## Script purpose: Global script for Shiny Tool
## Date: June 17th, 2017
## Author: Tyler Clavelle
##################################################

# Packages & Data ---------------------------------------------------------

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

### Source Functions
source('sardine_app_functions.R')
source('sardine_optim.R')

### Data
length_df <- read_csv(file = 'data/s_lemuru_length_freq.csv') # length frequency data
comm_df <- read_csv(file = 'data/C80PNVPPcm.csv', na = '...') # quarterly commercial landings data
muni_df <- read_csv(file = 'data/D10PNVMPcm.csv', na = '...') # quarterly municipal landings data
catch_annual <- read_csv(file = 'data/phils_fishery_compiled.csv') %>% # national catch data
  filter(CommName == 'Indian sardines (Tamban)')
# soi <- download_enso(climate_idx = "soi") # Southern Oscillation Index
soi <- read_csv(file = 'data/soi_timeseries.csv') # Southern Oscillation Index

vbK <- 1.29 # von Bert K
Linf <- 22.1 # L infinity
lwA <- 0.0078 # alpha from length-weight relationship
lwB <- 3.165 # beta from length weight relationship
t0 <- -lwA/lwB
M <- 1.786
price <- catch_annual$Price[nrow(catch_annual)]

length_at_age <- Linf * (1-exp(-(vbK/12)*(c(1:48) - t0)))
weight_at_age <- lwA * length_at_age ^ lwB


# Catch data processing for dygraph ---------------------------------------
##################################################
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


# Data processing for initial population ----------------------------------
##################################################
# Add weight
length_df <- length_df %>%
  filter(area == 'combined') %>%
  dplyr::select(length, freq) %>%
  mutate(weight = lwA * length ^ lwB)

# Final years catch value
init_f_weight <- catch_annual$Catch[1]

# Calculate average sardine weight from the length-frequency data and the age-length relationship
avg_weight <- length_df %>%
  mutate(total_wt = freq * weight) %>%
  summarize(avg_weight = sum(total_wt, na.rm = T) / sum(freq, na.rm = T)) %>% {.$avg_weight}

# Calculate total individuals by dividing catch (first convert to grams) by average weight
init_f_num <- (init_f_weight * 1e6) / avg_weight

# Break apart catch history into months
catch_history <- catch_annual %>% {.$Catch}
catch_history <- rep(catch_history / 12, each = 12)


# SOI Processing ----------------------------------------------------------

# Calculate annual SOI for optimization
soi_monthly <- soi %>%
  filter(Year %in% unique(catch_annual$Year))
