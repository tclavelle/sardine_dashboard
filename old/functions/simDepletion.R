#######################################################################################
## Project: EDF Sardine Dashboard
## Purpose: Calculate depletion for simulation
## Date: 11/15/2017
#######################################################################################

simDepletion <- function(mature_df, sim_name, weight_age, recruit_months, sim_length, ssb0_wt) {
  
  # Pull out recruitment months to calculate SSB
  mature_df <- mature_df * weight_age
  
  # Calculate depletion for each year based on the spawning stock biomass in the month of recruitment
  depletion_df <- mature_df %>%
    tbl_df() %>%
    mutate(month     = c(1:nrow(.)),
           year      = ceiling(rep_len(x = c(1:sim_length / 12), length.out = sim_length + 1)),
           total     = rowSums(.),
           depletion = 1 - total / ssb0_wt,
           scenario  = sim_name) %>%
    select(scenario, month, year, total, depletion) %>%
    filter(month %in% recruit_months)
  
  return(depletion_df)
}