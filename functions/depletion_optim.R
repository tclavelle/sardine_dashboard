#######################################################################################
## Project: EDF Sardine Dashboard
## Purpose: Fishery optimization function
## Date: 2017-10-16
#######################################################################################

depletion_NLL	<- function(par_start = c(1e6, 1), depletion, catch, soi) {
  
  # Simulate population from r0 out to equilibrium
  base_sim <- sardine_sim(r0 = par_start[1], 
                          f=0, 
                          recruit_month = 1, 
                          recruit_type = 'constant')
  
  # Set B0 to sum of spawning stock biomass in numbers
  B0 <- sum(base_sim[[6]], na.rm = T)
  ssb0 <- B0
  
  # Arbitrary value of R0
  R0 <-  par_start[1]
  
  # Assume steepness of 0.8
  h  <- 0.8

  # Bev Holt parameters
  alpha_bh		<-	(B0 / R0) * (1-h)/(4*h)
  beta_bh		<-	(5*h-1)/(4*h*R0)
  
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
