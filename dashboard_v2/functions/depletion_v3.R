#######################################################################################
## Project: EDF Sardine Dashboard
## Purpose: Optimize
## Date:
#######################################################################################

depletion_v3 <- function(par_start, target_depletion, r0 = 1e6, ssb0_wt, initial_pop, 
                         sim_length, weight_at_age, length_at_age, maturity,
                         recruit_months) {
  
  # Try forward simulation given F parameter 
  simForward <- sardine_sim_v3(initial_pop = initial_pop,
                               recruit_type = 'bev_holt',
                               r0 = r0,
                               ssb0 = ssb0_wt,
                               f = par_start,
                               ages = c(7:max_age),
                               sim_length = sim_length,
                               length_at_age = length_at_age,
                               weight_at_age = weight_at_age,
                               maturity = maturity)
  
  # Calculate depletion of simulation
  sim_dp <- simDepletion(mature_df = simForward[['m_out']], sim_name = 'current F', weight_age = weight_at_age, 
                             sim_length = sim_length, ssb0_wt = ssb0_wt, recruit_months = recruit_months)
  
  # calculate sum of squares value to minimize
  out <- sim_dp %>%
    mutate(diff = (depletion - target_depletion)^2) %>%
    summarize(sum_of_squares = sum(diff, na.rm = T)) %>% .$sum_of_squares
  
  return(out)

}