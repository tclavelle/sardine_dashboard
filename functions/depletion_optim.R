#######################################################################################
## Project: EDF Sardine Dashboard
## Purpose: Fishery optimization function
## Date: 2017-10-16
#######################################################################################

depletion_NLL	<- function(
  par_start = c(1e6, 1),
  depletion, 
  catch, 
  soi) {
  # browser()
  base_sim <- sardine_sim(r0 = par_start[1], 
                          f=0, 
                          recruit_month = 1, 
                          recruit_type = 'constant')
  
  # Set B0 to sum of spawning stock biomass in numbers
  B0 <- sum(base_sim[[6]], na.rm = T)
  # Pull out SSB in MT
  # ssb0 <- base_sim[[5]]
  ssb0 <- B0
  
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
  # browser()
  test <- fishing_sim$f_timeseries
  out_depletion <- fishing_sim$final_depletion
  
  return(sum((out_depletion - depletion)^2))
}
