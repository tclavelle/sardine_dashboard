# ssb0
# m = 1.786
# vbk = 1.29
# t0 = -0.08
# lwA = 0.0078
# lwB = 3.165
# linf = 22.1
# r0 = 1e10
# selectivity = 'selA'
# catch
# ages = c(7:48) # month 7 is age of recruitment to fishery
# age_mat = c(7:10) # months over which fish attain maturity
# recruit_age = 7
# recruit_month = 1
# sim_length = 156
# sex_ratio = 0.5
# recruit_type = 'bev_holt'
# alpha_bh
# beta_bh
# closed <- NA

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
                        recruit_var = 1,
                        soi_project = 1,
                        alpha_bh,
                        beta_bh,
                        closed = NA) {
  
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
  # browser()
  ## Set harvest regime for whole timeseries
  adult_harvest   <- rep(1- exp(-1/12 * f),nrow(n_out))
  # Set F in closed season to 0
  adult_harvest[closed] <- 0
  
  # Set percent survivors. This value represents total natural mortality per month. 
  p_surv <- exp(- (1 / 12 * (m)))
  p_surv	<- rep(p_surv,nrow(n_out)) 
  
  ### Recruitment settings
  # Function for turning CV on the log-normal scale into standard deviations you can plug into R functions 
  log_norm_cv		<- function(cv){
    sigma2	<-	 log(cv^2 + 1)
    return(sqrt(sigma2))
  }
  
  # Set recruitment months
  recruit_months <- seq(from = recruit_month, to = sim_length, by = 12)
  
  # Section: Simulate fishery through time in monthly timesteps
  #######################################################################################				
  for(i in 1:(sim_length)){
    ## Recruitment
    
    # constant recruitment
    if(recruit_type == 'constant') { 
      if(i %in% recruit_months) {
        n_out[i+1,recruit_age] <- r0 
      } else n_out[i+1,recruit_age] <- 0 
    }
    # Beverton-Holt 
    if(recruit_type == 'bev_holt') {
      # browser()
      # if(i == 19) {browser()}
      if(i %in% recruit_months) {
        
        n_out[i+1,recruit_age]	<-	sex_ratio * bev_holt(alpha_bh,
                                                       beta_bh,
                                                       sum(m_out[i,recruit_age:ncol(m_out)], na.rm = T)) * recruit_var[i] * soi_project[i]
      } else n_out[i+1,recruit_age] <- 0
    }
    # browser()
    # Survive to next year	
    n_out[(i+1),(recruit_age + 1):ncol(n_out)]	<- 	n_out[i, recruit_age:c(ncol(n_out)-1)] * p_surv[i]
    n_out[(i+1),ncol(n_out)]		<-	n_out[i,(ncol(n_out)-1)] * p_surv[i] + n_out[i,ncol(n_out)] * p_surv[i]
    
    # Calculate Biomass
    b_out[i+1,]	<-	(weight_at_age * n_out[i+1,]) / 1e6 # divide by 1e6 to convert grams to metric tons
    
    # Mature individuals 
    m_out[i+1,] <- n_out[i+1,] * maturity
    
    # If there's a fishery
    if(adult_harvest[i+1] > 0){
      # browser()
      # Determine Harvest Rate
      tot_harvest 	<-	 adult_harvest[i+1] * sum(m_out[i+1,])
       
      # Fishery Happens on Mature Fish
      effect_harvest_rate		<- tot_harvest/sum(m_out[i+1,] * select)
      
      if(f_mode == 'rate') { c_out[i+1,] <- effect_harvest_rate * m_out[i+1,] * select }
      if(f_mode == 'catch') { c_out[i+1,] <- f[i] * select * weight_at_age }
      
      # Calibrate Mature individuals
      m_out[i+1,] <- m_out[i+1,] - c_out[i+1,]
      
      # Calibrate total biomass
      b_out[i+1,]		<- b_out[i+1,] - c_out[i+1,] * weight_at_age / 1e6
      
      # Calibrate N_individuals
      # n_out[i+1,]		<-	b_out[i+1,] * 1e6 / weight_at_age
    }
  } # close year loop

  # Initial pawning stock biomass (K)
  K <- sum(m_out[recruit_months[length(recruit_months)],ages] * weight_at_age / 1e6, na.rm = T)
  K_num <- m_out[recruit_months[length(recruit_months)],ages]
  
  # Section: Add labels to result data frames
  #######################################################################################
  n_out <- n_out %>%
    tbl_df() %>%
    mutate(month = c(1:(sim_length+1)),
           year  = ceiling(rep_len(x = c(1:sim_length /12), length.out = sim_length + 1))) %>%
    select(month, year, everything()) %>%
    gather(key = 'age_class', value = 'abundance', 3:ncol(.))
  
  b_out <- b_out %>%
    tbl_df() %>%
    mutate(month = c(1:(sim_length+1)),
           year  = ceiling(rep_len(x = c(1:sim_length /12), length.out = sim_length + 1))) %>%
    select(month, year, everything()) %>%
    gather(key = 'age_class', value = 'biomass', 3:ncol(.))
  
  m_out <- m_out %>%
    tbl_df() %>%
    mutate(month = c(1:(sim_length+1)),
           year  = ceiling(rep_len(x = c(1:sim_length /12), length.out = sim_length + 1))) %>%
    select(month, year, everything()) %>%
    gather(key = 'age_class', value = 'mature', 3:ncol(.))
  
  # browser()
  c_out <- c_out * weight_at_age / 1e6
  c_out <- c_out %>%
    tbl_df() %>%
    mutate(month = c(1:(sim_length+1)),
           year  = ceiling(rep_len(x = c(1:sim_length /12), length.out = sim_length + 1))) %>%
    select(month, year, everything()) %>%
    gather(key = 'age_class', value = 'catch', 3:ncol(.)) %>%
    mutate(revenue = catch * price)
  
  return(list(n_out, b_out, c_out, m_out, K, K_num))
} 


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
  # browser()
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
  n_out[1,1:length(r0)] <- r0
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
    b_out[i+1,]	<-	(weight_at_age * n_out[i+1,]) / 1e6 # divide by 1e6 to convert grams to metric tons
    
    # Mature individuals 
    m_out[i+1,] <- n_out[i+1,] * maturity * 0.5
    
    # If there's a fishery
    if(catch[i] > 0){
      
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
 
  return(list(abundance = n_out, biomass = b_out, catch = c_out, final_depletion = final_depletion, ssb = b_final))
}


# Optimization ---------------------------------------------------

# Minimization function
depletion_NLL	<- function(
  par_start = c(1e6, 1),
  depletion, 
  catch, 
  soi){
  # browser()
  base_sim <- sardine_sim(r0 = par_start[1], 
                          f=0, 
                          recruit_month = 1, 
                          recruit_type = 'constant')
  
  # Pull out SSB in MT
  ssb0 <- base_sim[[5]]
  
  # Set B0 to sum of spawning stock biomass in numbers
  B0 <- sum(base_sim[[6]], na.rm = T)
  # Arbitrary value of R0
  R0 <-  par_start[1]
  # Assume steepness of 0.8
  h  <- 0.8
  # Temporary values to optimize
  alpha.temp		<-	(B0 / R0) * (1-h)/(4*h)
  beta.temp		<-	(5*h-1)/(4*h*R0)		
  
  R	<-	bev_holt(alpha.temp,beta.temp,(seq(1, B0, by=B0/100)))
  
  alpha_bh <- alpha.temp
  beta_bh <- beta.temp
  
  # Adjust SOI by SOI parameter value
  soi_for_opt <- soi * par_start[2]
  
  fishing_sim <- sardine_optim(ssb0 = ssb0,
                               r0 = R0,
                               catch = catch, 
                               soi = soi_for_opt, 
                               recruit_month = 1, 
                               alpha_bh = alpha_bh, 
                               beta_bh = beta_bh, 
                               recruit_type = 'bev_holt')
  
  out_depletion <- fishing_sim$final_depletion
  
  return(sum((out_depletion - depletion)^2))
}
