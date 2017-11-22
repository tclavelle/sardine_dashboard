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
    # geom_smooth(method = 'loess') +
    scale_y_continuous(labels = comma) +
    labs(x = 'Time (months)',
         y = metric) +
    theme_bw()
    
  return(plot_out)
}

# Simulation depletion plot
depletion_plot <- function(depletion_df) {
  
  dep_plot_out <- depletion_df %>%
    ggplot(aes(x = year, y = depletion, color = scenario)) +
    geom_line() +
    # geom_smooth(se = F) +
    labs(x = 'Time (year)',
         y = 'Depletion') +
    theme_bw()
  
  return(dep_plot_out)
}

# Simulation timeseries comparison plot
simCompare <- function(simA, simB, sim_length, metric) {
  
  # Tidy both data frames
  simA <- simA %>%
    tbl_df() %>%
    mutate(month = c(1:nrow(.)),
           year = ceiling(rep_len(x = c(1:sim_length / 12), length.out = sim_length + 1))) %>%
    gather(key = 'age', value = 'number', 1:c(ncol(.)-2)) %>%
    group_by(year, month) %>%
    summarize(statusQuo = sum(number, na.rm = T)) 
  
  simB <- simB %>%
    tbl_df() %>%
    mutate(month = c(1:nrow(.)),
           year = ceiling(rep_len(x = c(1:sim_length / 12), length.out = sim_length + 1))) %>%
    gather(key = 'age', value = 'number', 1:c(ncol(.)-2)) %>%
    group_by(year, month) %>%
    summarize(intervention = sum(number, na.rm = T)) 
  
  # Bind data frames
  compare_df <- left_join(simA, simB) %>%
    mutate(difference = 100 * (intervention - statusQuo) / statusQuo)
  
  # Plot
  compare_out <- compare_df %>%
    ggplot(aes(x = month, y = difference)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 2, color = 'red') +
    # geom_smooth(method = 'loess') +
    scale_y_continuous(labels = comma) +
    labs(x = 'Time (months)',
         y = metric,
         title = paste0(metric, ' relative to status quo')) +
    theme_bw()
  
  return(compare_out)
}

