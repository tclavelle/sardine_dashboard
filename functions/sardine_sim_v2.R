#######################################################################################
## Project: EDF Sardine Dashboard
## Purpose: Fishery simulation function
## Date: 2017-10-16
#######################################################################################

sardine_sim_v2 <- function(ssb0 = 6188691,
                         r0 = 1e5, 
                         m = 1.786,
                         f = 0,
                         length_at_age,
                         weight_at_age,
                         maturity,
                         selectivity = 'selA', 
                         f_mode = 'catch',
                         ages = c(1:42), # month 7 is age of recruitment to fishery
                         age_mat = c(1:4), # months over which fish attain maturity
                         recruit_month = 1,
                         sim_length = 100,
                         recruit_type = 'constant',
                         recruit_var = 1,
                         closed = NA) {
  

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
  
  colnames(n_out) <- c(paste0(ages + 6, '_month_olds')) # set column names
  b_out <- n_out # biomass data frame
  m_out <- n_out # mature biomass data frame
  c_out <- n_out # catch data frame
  
  # F timeseries data frame
  f_out <- data_frame(month = c(1:sim_length),
                      f     = rep.int(0, times = sim_length))
  
  ## Set initial conditions
  n_out[1,1:length(r0)] <- r0
  b_out[1,] <- n_out[1,] * weight_at_age / 1e6 # divide by 1e6 to put convert biomass from grams to metric tons
  m_out[1,] <- n_out[1,] * maturity # multiply number of individual times percent mature
  
  ## Set harvest regime for whole timeseries
  adult_harvest   <- rep(1 - exp(-(1/12 * f)),nrow(n_out))
  # Set F in closed season to 0
  adult_harvest[closed] <- 0
  
  # Set percent survivors. This value represents total natural mortality per month. 
  p_surv <- exp(- (1 / 12 * (m)))
  p_surv	<- rep(p_surv,nrow(n_out)) 
  
  # Set recruitment months
  recruit_months <- seq(from = recruit_month, to = sim_length, by = 12)
  # browser()
  # Section: Simulate fishery through time in monthly timesteps
  #######################################################################################				
  for(i in 1:(sim_length)){
    
    # Constant recruitment
    if(recruit_type == 'constant') { 
      if(i %in% recruit_months) {
        n_out[i+1,1] <- r0 # constant recruitment of R0 in recruitment months
      } else n_out[i+1,1] <- 0 # otherwise zero recruitment
    }
    
    # Beverton-Holt 
    if(recruit_type == 'bev_holt') {
      
      if(i %in% recruit_months) {
        
        ssb <- sum(m_out[i,] * weight_at_age)
        n_out[i+1,1]	<- (0.8 * r0 * 0.8 * ssb) / (0.2 * ssb0 * (1 - 0.8) + (0.8 - 0.2) * ssb)	
          
      } else n_out[i+1,1] <- 0
    }
    
    # Survive to next year	
    n_out[(i+1),2:ncol(n_out)]	<- 	n_out[i, 1:c(ncol(n_out)-1)] * p_surv[i]
    n_out[(i+1),ncol(n_out)]		<-	n_out[i,(ncol(n_out)-1)] * p_surv[i] + n_out[i,ncol(n_out)] * p_surv[i]
    
    # Calculate Biomass
    b_out[i+1,]	<-	weight_at_age * n_out[i+1,] / 1e6 # divide by 1e6 to convert grams to metric tons
    
    # Mature individuals 
    m_out[i+1,] <- n_out[i+1,] * maturity #* sex_ratio
    
    # If there's a fishery
    if(adult_harvest[i+1] > 0){
      
      # Determine Harvest Rate
      tot_harvest 	<-	 adult_harvest[i+1] * sum(m_out[i+1,])
      # Fishery Happens on Mature Fish
      effect_harvest_rate		<- tot_harvest/sum(m_out[i+1,] * select)
      
      if(f_mode == 'rate') { c_out[i+1,] <- effect_harvest_rate * m_out[i+1,] * select }
      
      if(f_mode == 'catch') { 
        # Use catch history data
        c_out[i+1,] <- f[i] * select * weight_at_age 
        # Calculate F and save timeseries
        f_out[i,2] <- f[i] / b_out[i]
      }
      
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
  
  return(list(n_out = n_out, b_out = b_out, c_out = c_out, m_out = m_out, K = K, K_num = K_num, f_out = f_out))
} 