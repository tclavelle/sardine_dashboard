#######################################################################################
## Project: EDF Philippines Sardine Simulation Functions
## Purpose: Functions for Sardine web simulation app
## Date: April 16th, 2017
## Author: Tyler Clavelle
#######################################################################################

# Probability of selection function
gill_sel <- function(l,lopt){
  exp(-(l-lopt)^2/(2*1.747^2))
  }

# r0 = 4e10
# m = 1.786
# vbk = 1.29
# t0 = -0.08
# lwA = 0.0078
# lwB = 3.165
# linf = 22.1
# f = 0.5
# init_f_number = c(4e10)
# selectivity = 'selA'
# ages = c(1:48)
# age_mat = c(7:10)
# recruit_age = 7
# sim_length = 100
# sex_ratio = 0.5
# recruit_type = 'bev_holt'
# alpha_bh
# beta_bh

## !! Write two separate functions - 1 to calibrate the model, and one for the simulation

sardine_sim <- function(r0 = 4e10, 
                        m = 1.786,
                        vbk = 1.29,
                        t0 = -0.08,
                        lwA = 0.0078,
                        lwB = 3.165,
                        linf = 22.1,
                        f = 0,
                        init_f_number = 4e10,
                        selectivity = 'selA', 
                        f_mode = 'catch',
                        ages = c(1:48), # month 7 is age of recruitment to fishery
                        age_mat = c(7:10), # months over which fish attain maturity
                        recruit_age = 7,
                        recruit_month = 1,
                        sim_length = 100,
                        sex_ratio = 0.5,
                        recruit_type = 'constant',
                        alpha_bh,
                        beta_bh) {
  
  # generate length and weight at age vectors for number of ages
  length_at_age <- linf * (1-exp(-(vbk/12)*(ages - t0)))
  weight_at_age <- lwA * length_at_age ^ lwB
  # browser()
  # generate maturity vector
  mat_ramp <- seq(0,1, length.out = length(age_mat)+1) # assume proportion mature is linear between ages
  maturity <- c(rep(0, times = age_mat[1]-1), # zero mature before length at first mature (11.5 cm) 
                mat_ramp[2:length(mat_ramp)], # maturity ramp
                rep(1, times = length(length_at_age)-max(age_mat))) # all large fish are mature
  
  # Probability of selection function
  gill_sel <- function(l,lopt){exp(-(l-lopt)^2/(2*1.747^2))}
  # Selectivity of 2.54cm mesh
  selA <- sapply(length_at_age, gill_sel, lopt = 14.96)
  # Selectivity of 3.175cm mesh
  selB <- sapply(length_at_age, gill_sel, lopt = 18.7)
  
  if(selectivity == 'selA') { select = selA } else select = selB
  
  # Section: Initialize data frames for storing simulation results
  #######################################################################################
  ## Build data frames
  n_out <- matrix(0,ncol = length(ages), # extra column to track months in simulation
                             nrow = sim_length + 1) 
  colnames(n_out) <- c(paste0(ages, '_month_olds')) # set column names
  b_out <- n_out # biomass data frame
  m_out <- n_out # mature biomass data frame
  c_out <- n_out # catch data frame
  
  ## Set initial conditions
  n_out[1,recruit_age:(recruit_age + (length(init_f_number)-1))] <- init_f_number
  b_out[1,] <- n_out[1,] * weight_at_age / 1e6 # divide by 1e6 to put convert biomass from grams to metric tons
  m_out[1,] <- n_out[1,] * maturity * sex_ratio # multiply number of individual times percent mature
  
  ## Set harvest regime for whole timeseries
  adult_harvest   <- rep(1- exp(-1/12 * f),nrow(n_out))
  
  # Set percent survivors. This value represents total natural mortality per month. 
  p_surv <- exp(- (1 / 12 * (m + f)))
  p_surv	<- rep(p_surv,nrow(n_out)) 
  
  ## Set recruitment months
  recruit_months <- seq(from = recruit_month, to = sim_length, by = 12)
  
  # Section: Simulate fishery through time in monthly timesteps
  #######################################################################################				
  for(i in 1:(sim_length)){
    # print(i)
    ## Recruitment
    # constant recruitment
    if(recruit_type == 'constant') { 
      if(i %in% recruit_months) {
        n_out[i+1,recruit_age] <- r0 
      } else n_out[i+1,recruit_age] <- 0 
    }
    # Beverton-Holt 
    if(recruit_type == 'bev_holt') {
      if(i %in% recruit_months) {
        n_out[i+1,recruit_age]	<-	sex_ratio * bev_holt(alpha_bh,beta_bh,sum(m_out[i,recruit_age:ncol(m_out)], na.rm = T))
    } else n_out[i+1,recruit_age] <- 0
    }
    
    # Survive to next year	
    n_out[(i+1),(recruit_age + 1):ncol(n_out)]	<- 	n_out[ i,(recruit_age):c(ncol(n_out)-1)]*p_surv[i]
    n_out[(i+1),ncol(n_out)]		<-	n_out[ i,(ncol(n_out)-1)]*p_surv[i] + n_out[ i,ncol(n_out)]*p_surv[i]
    
    # Calculate Biomass
    b_out[i+1,]	<-	(weight_at_age * n_out[i+1,]) / 1e6 # divide by 1e6 to convert grams to metric tons
    
    # Mature individuals 
    m_out[i+1,] <- n_out[i+1,] * maturity * 0.5
    
    # If there's a fishery
    if(adult_harvest[i] > 0){
      # Determine Harvest Rate
      tot_harvest 	<-	 adult_harvest[i+1] * sum(m_out[i+1,])
      # 
      # # Fishery Happens on Mature Fish
      effect_harvest_rate		<- tot_harvest/sum(m_out[i+1,] * select)
      # c_out[i+1,]			<- effect_harvest_rate * m_out[i+1,] * select
      if(f_mode == 'rate') { c_out[i+1,] <- effect_harvest_rate * m_out[i+1,] * select }
      if(f_mode == 'catch') { c_out[i+1,] <- f[i] * select * weight_at_age }
      
      # Calibrate Mature individuals
      m_out[i+1,] <- m_out[i+1,] - c_out[i+1,]
      
      # Calibrate total biomass
      b_out[i+1,]		<- b_out[i+1,] - (c_out[i+1,] * weight_at_age / 1e6)
      
      # Calibrate N_individuals
      n_out[i+1,]		<-	b_out[i+1,] * 1e6 / weight_at_age
    }
  } # close year loop
  
  # Initial pawning stock biomass (K)
  K <- sum(m_out[recruit_months[length(recruit_months)],ages] * weight_at_age / 1e6, na.rm = T)
  K_num <- m_out[recruit_months[length(recruit_months)],ages]
  
  # Section: Add labels to result data frames
  #######################################################################################
  n_out <- n_out %>%
    tbl_df() %>%
    mutate(month = c(1:(sim_length+1))) %>%
    select(month, everything()) %>%
    gather(key = 'age_class', value = 'abundance', 2:ncol(.))
  
  b_out <- b_out %>%
    tbl_df() %>%
    mutate(month = c(1:(sim_length+1))) %>%
    select(month, everything()) %>%
    gather(key = 'age_class', value = 'biomass', 2:ncol(.))
  
  m_out <- m_out %>%
    tbl_df() %>%
    mutate(month = c(1:(sim_length+1))) %>%
    select(month, everything()) %>%
    gather(key = 'age_class', value = 'mature', 2:ncol(.))
  
  c_out <- c_out %>%
    tbl_df() %>%
    mutate(month = c(1:(sim_length+1))) %>%
    select(month, everything()) %>%
    gather(key = 'age_class', value = 'catch', 2:ncol(.)) 
  
  return(list(n_out, b_out, c_out, m_out, K, K_num))
} 


# Probability of selection function
gill_sel <- function(l,lopt){exp(-(l-lopt)^2/(2*1.747^2))}

# Beverton-Holt
bev_holt	<- 	function(alpha,beta,b){
  rec <- b * (alpha + beta * b) ^ (-1)
  return(rec)
}

# Negative log-likelihood fucntion
NLL	<- function(PAR.start,Rec, B){
  ALPHA = PAR.start[1]
  BETA  = PAR.start[2]
  Recruits = bev_holt(ALPHA,BETA,B)
  
  return(sum((Recruits - Rec)^2))
}

# Logit model for estimating probability of capture (selectivity)
Logit	<- function(a,b,x){	1 / (1+exp(-(a+b*x)))}

# Function for turning CV on the log-normal scale into standard deviations you can plug into R functions 
log.norm.cv		<- function(cv){
  sigma2	<-	 log(cv^2 + 1)
  return(sqrt(sigma2))
}

