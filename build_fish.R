#######################################################################################
## Project: EDF Sardine Dashboard
## Purpose: Fishery simulation function
## Date: 2017-10-16
#######################################################################################

create_fish <- function(common_name = 'Bali sardine',
                        scientific_name = "Sardinella lemuru",
                        linf = NA,
                        vbk = NA,
                        t0 = -0.08,
                        length_units = 'cm',
                        max_age = 20,
                        weight_a = 0.0078,
                        weight_b = 3.165,
                        weight_units = 'kg',
                        length_50_mature = NA,
                        length_95_mature = NA,
                        age_50_mature = NA,
                        age_95_mature = NA,
                        age_mature = NA,
                        length_mature = 12.5,
                        m = 1.786,
                        steepness = 0.8,
                        r0 = 1000,
                        price = 1) {
  
  fish <- list()

  fish$length_at_age <- linf * (1 - exp(-vbk * ((1:max_age - t0))))
  
  # process weight
  fish$weight_at_age <- weight_a * fish$length_at_age ^ weight_b
  
  # process maturity
  if ((is.na(age_50_mature) |
       is.na(age_95_mature)) & is.na(age_mature) == F) {
    
    age_50_mature <- age_mature
    age_95_mature <- 1.01 * age_50_mature
    
  } else if (is.na(age_mature)) {
    if (is.na(length_mature)) {
      length_mature <-  linf * lmat_to_linf_ratio
    }
    
    age_mature <- (log(1 - length_mature / linf) / -vbk) + t0
    
    age_50_mature <- age_mature
    
    age_95_mature <- 1.01 * age_50_mature
  }
  
  fish$maturity_at_age <-
    ((1 / (1 + exp(-log(
      19
    ) * (((1:max_age) - age_50_mature) / (age_95_mature - age_50_mature)
    )))))
  
  fish$scientific_name <- scientific_name
  fish$common_name <- common_name
  fish$ssb_at_age <- fish$maturity_at_age * fish$weight_at_age
  fish$linf <- linf
  fish$vbk  <-  vbk
  fish$t0 <-  t0
  fish$max_age <-  max_age
  fish$weight_a <-  weight_a
  fish$weight_b <-  weight_b
  fish$length_50_mature <-  length_50_mature
  fish$length_95_mature <-  length_95_mature
  fish$age_50_mature <-  age_50_mature
  fish$age_95_mature <-  age_95_mature
  fish$age_mature <-  age_mature
  fish$length_mature <-  length_mature
  fish$m <-  m
  fish$steepness <- steepness
  fish$r0 <- r0
  fish$density_dependence_form = density_dependence_form
  fish$adult_movement <-  adult_movement
  fish$larval_movement <-  larval_movement
  fish$lmat_to_linf_ratio <-  lmat_to_linf_ratio
  fish$length_units <-  length_units
  fish$weight_units <-  weight_units
  fish$price <- price
  
  return(fish)
}