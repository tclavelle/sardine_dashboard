#######################################################################################
## Project: EDF Sardine Dashboard
## Purpose: Fishery simulation function
## Date: 2017-10-16
#######################################################################################

sardine_optim_v2 <- function(r0 = 1e5,
                             ssb0 = 6188691, 
                             m = 1.786,
                             length_at_age,
                             weight_at_age,
                             maturity,
                             selectivity = 'selA', 
                             catch,
                             ages = c(7:48), # month 7 is age of recruitment to fishery
                             age_mat = c(1:4), # months over which fish attain maturity
                             recruit_month = 1,
                             sim_length = 156,
                             recruit_type = 'bev_holt') {
  
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
  
  colnames(n_out) <- c(paste0(ages, '_month_olds')) 
  b_out <- n_out # biomass data frame
  m_out <- n_out # mature biomass data frame
  c_out <- n_out # catch data frame
  
  ## Set initial conditions
  n_out[1,1:length(r0)] <- r0
  b_out[1,] <- n_out[1,] * weight_at_age / 1e6 # divide by 1e6 to put convert biomass from grams to metric tons
  m_out[1,] <- n_out[1,] * maturity  # multiply number of individual times percent mature
  
  ## Set recruitment months
  recruit_months <- seq(from = recruit_month, to = sim_length, by = 12)
  
  # Section: Simulate fishery through time in monthly timesteps
  #######################################################################################				
  for(i in 1:(sim_length)){
    
    ## Recruitment
    # constant recruitment
    if(recruit_type == 'constant') { 
      if(i %in% recruit_months) {
        n_out[i+1,1] <- r0 
      } else n_out[i+1,1] <- 0
    }
    
    # Beverton-Holt 
    if(recruit_type == 'bev_holt') {
      
      if(i %in% recruit_months) {
        # browser()
        ssb <- sum(m_out[i,] * weight_at_age)
        n_out[i+1,1]	<- (0.8 * r0 * 0.8 * ssb) / (0.2 * ssb0 * (1 - 0.8) + (0.8 - 0.2) * ssb)	
        
      } else n_out[i+1,1] <- 0
    }
    
    # Calculate highest possible F given available biomass
    # highest_f <- sum(n_out[i+1,] * select * weight_at_age) / sum(m_out[i,] * weight_at_age)
    # catch_f <- catch[i] / (sum(m_out[i,] * weight_at_age) / 1e6)
    # 
    # if(catch_f > highest_f) { 
    #   f <- highest_f
    # } else f <- catch_f
    
    # Set survival
    p_surv <- exp(- (1 / 12 * (m)))
    
    # assign catch
    # c_out[i+1,] <- n_out[i+1,] * select * f
    
    # Survive to next year	
    n_out[(i+1),2:ncol(n_out)]	<- 	n_out[i,1:c(ncol(n_out)-1)]*p_surv
    n_out[(i+1),ncol(n_out)]		<-	n_out[ i,(ncol(n_out)-1)]*p_surv + n_out[ i,ncol(n_out)]*p_surv
    
    # Calculate Biomass
    b_out[i+1,]	<-	weight_at_age * n_out[i+1,] / 1e6 # divide by 1e6 to convert grams to metric tons
    
    # Mature individuals 
    m_out[i+1,] <- n_out[i+1,] * maturity 
    
    # If there's a fishery
    if(catch[i] > 0){

      # Calculate F from catch on mature fish
      f <- catch[i] / sum(m_out[i+1,] * weight_at_age / 1e6)

      # What the catch would be given the calculated F
      temp_catch <- f * n_out[i+1,] * select

      # Where the catch would be larger than available biomass, set to available biomass
      # temp_catch[temp_catch > m_out[i+1,]] <- m_out[i+1,][temp_catch > m_out[i+1,]]

      # Fishery Happens on Mature Fish
      c_out[i+1,] <- temp_catch

      # Calibrate Mature individuals
      m_out[i+1,] <- m_out[i+1,] - c_out[i+1,]

      # Calibrate total biomass
      b_out[i+1,]		<- b_out[i+1,] - (c_out[i+1,] * weight_at_age / 1e6)

      # Calibrate N_individuals
      n_out[i+1,]		<-	b_out[i+1,] * 1e6 / weight_at_age
    }
  } # close year loop
  
  # Section: Add labels to result data frames
  #######################################################################################
  
  b_final <- sum(m_out[nrow(m_out),] * weight_at_age, na.rm = T)
  
  # Calculate final depletion
  final_depletion <- b_final / ssb0 
  
  return(list(abundance = n_out, 
              biomass = b_out, 
              catch = c_out, 
              final_depletion = final_depletion, 
              ssb = b_final,
              f_timeseries = data_frame(year = c(1:sim_length),
                                        f    = f)))
}
