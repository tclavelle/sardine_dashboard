#######################################################################################
## Project: EDF Sardine Dashboard
## Purpose: Fishery simulation function
## Date: 2017-10-16
#######################################################################################

sardine_optim <- function(r0 = 1e10,
                          ssb0, 
                          m = 1.786,
                          vbk = 1.29,
                          t0 = -0.08,
                          lwA = 0.0078,
                          lwB = 3.165,
                          linf = 22.1,
                          selectivity = 'selA', 
                          catch,
                          soi,
                          ages = c(7:48), # month 7 is age of recruitment to fishery
                          age_mat = c(7:10), # months over which fish attain maturity
                          recruit_month = 1,
                          sim_length = 156,
                          sex_ratio = 0.5,
                          recruit_type = 'bev_holt',
                          alpha_bh,
                          beta_bh) {
  
  # generate length and weight at age vectors for number of ages
  length_at_age <- linf * (1-exp(-(vbk/12)*(ages - t0)))
  weight_at_age <- lwA * length_at_age ^ lwB
  
  # generate maturity vector
  mat_ramp <- seq(0,1, length.out = length(age_mat)+1) # assume proportion mature is linear between ages
  maturity <- c(mat_ramp[2:length(mat_ramp)], # maturity ramp
                rep(1, times = length(length_at_age)-(length(mat_ramp)-1))) # all large fish are mature
  
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
  n_out[1,1] <- r0
  b_out[1,] <- n_out[1,] * weight_at_age / 1e6 # divide by 1e6 to put convert biomass from grams to metric tons
  m_out[1,] <- n_out[1,] * maturity * sex_ratio # multiply number of individual times percent mature
  
  ## Set recruitment months
  recruit_months <- seq(from = recruit_month, to = sim_length, by = 12)
  
  # Section: Simulate fishery through time in monthly timesteps
  #######################################################################################				
  for(i in 1:(sim_length)){
    
    ## Recruitment
    # constant recruitment
    if(recruit_type == 'constant') { n_out[i+1,1] <- r0 }
    # Beverton-Holt 
    if(recruit_type == 'bev_holt') {
      
      if(i %in% recruit_months) {
        # calculate number of recruits
        recruits <- sex_ratio * bev_holt(alpha_bh,beta_bh,sum(m_out[i,1], na.rm = T)) * 1e6
        # multiply times soi parameter
        n_out[i+1,1]	<-	recruits * soi[i]
      } else n_out[i+1,1] <- 0
      
    }
    
    # Set survival
    p_surv <- exp(- (1 / 12 * (m)))
    
    # Survive to next year	
    n_out[(i+1),2:ncol(n_out)]	<- 	n_out[i,1:c(ncol(n_out)-1)]*p_surv
    n_out[(i+1),ncol(n_out)]		<-	n_out[ i,(ncol(n_out)-1)]*p_surv + n_out[ i,ncol(n_out)]*p_surv
    
    # Calculate Biomass
    b_out[i+1,]	<-	weight_at_age * n_out[i+1,] / 1e6 # divide by 1e6 to convert grams to metric tons
    
    # Mature individuals 
    m_out[i+1,] <- n_out[i+1,] * maturity * sex_ratio
    
    # If there's a fishery
    if(catch[i] > 0){
      # browser()
      # Calculate F from catch on mature fish
      f <- catch[i] / sum(m_out[i+1,] * weight_at_age / 1e6) 
      
      # # Fishery Happens on Mature Fish
      c_out[i+1,] <- f * m_out[i+1,] * select 
      
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
 
  b_final <- sum(m_out[nrow(m_out),], na.rm = T)
 
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
