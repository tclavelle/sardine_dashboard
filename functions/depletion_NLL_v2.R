#######################################################################################
## Project: EDF Sardine Dashboard
## Purpose: Fishery optimization function
## Date: 2017-10-16
#######################################################################################

depletion_NLL_v2	<- function(par_start, depletion, catch, length_at_age, weight_at_age, maturity) {
  
  # Simulate population from r0 out to equilibrium
  base_sim <- sardine_sim_v2(r0 = par_start[1],
                             length_at_age = length_at_age,
                             weight_at_age = weight_at_age,
                             maturity = maturity,
                             f=0, 
                             recruit_month = 1, 
                             recruit_type = 'constant')
  
  # Generate virgin biomass from a given R0 (making sure to set ages in between recruitment events to 0)
  virgin <- par_start[1] * exp(- (1 / 12) * c(0:41))
  
  # Calculate number mature and multiply times weight at age to get SSB0
  ssb0 <- virgin * maturity * weight_at_age
  
  recruit_months <- seq(from = 1, to = 42, by = 12)
  
  # Only include numbers in recruitment 
  ssb0 <- sum(ssb0[recruit_months]) 
  
  fishing_sim <- sardine_optim_v2(ssb0 = ssb0,
                               r0 = par_start[1],
                               length_at_age = length_at_age,
                               weight_at_age = weight_at_age,
                               maturity = maturity,
                               catch = catch_history, 
                               recruit_month = 1, 
                               recruit_type = 'bev_holt')
  
  out_depletion <- fishing_sim$final_depletion
  
  if(out_depletion < 0) out_depletion <- 999
  
  return((out_depletion - depletion)^2)
}
