#######################################################################################
## Project: Sardine Dashboard
## Purpose: Simulation Plotting Function
## Date: 11-06-2017
#######################################################################################

# Simulation timeseries plots
sim_plot <- function(sim_df, sim_length, metric) {
  
  plot_out <- sim_df %>%
    tbl_df() %>%
    mutate(month = c(1:nrow(.)),
           year = ceiling(rep_len(x = c(1:sim_length / 12), length.out = sim_length + 1))) %>%
    gather(key = 'age', value = 'number', 1:c(ncol(.)-2)) %>%
    group_by(year, month) %>%
    summarize(Value = sum(number, na.rm = T)) %>%
    ggplot(aes(x = month, y = Value)) +
    geom_line() +
    geom_smooth(method = 'loess') +
    scale_y_continuous(labels = comma) +
    labs(x = 'Time (months)',
         y = metric) +
    theme_bw()
    
  return(plot_out)
}

# Simulation depletion plot
depletion_plot <- function(mature_df, weight_at_age, sim_length, spawning_bio) {
  
  # convert maturity data frame to weight
  mature_df <- mature_df * weight_at_age
  
  dep_plot_out <- mature_df %>%
    tbl_df() %>%
    mutate(month = c(1:nrow(.)),
           year = ceiling(rep_len(x = c(1:sim_length / 12), length.out = sim_length + 1)),
           total = rowSums(.),
           depletion = 1 - total / spawning_bio) %>%
    select(month, year, total, depletion) %>%
    ggplot(aes(x = month, y = depletion)) +
    geom_line() +
    geom_smooth(method = 'loess') +
    labs(x = 'Time (months)',
         y = 'Depletion') +
    theme_bw()
  
  return(dep_plot_out)
}