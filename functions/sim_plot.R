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
    scale_color_brewer(type = 'qual', palette = 'Dark2') +
    labs(x = 'Time (year)',
         y = 'Depletion') +
    theme_bw()
  
  return(dep_plot_out)
}

# Simulation timeseries comparison plot
simCompare <- function(simA, simB, sim_length, metric, sample_size) {
  
  simA <- simA %>%
    bind_rows() %>%
    rename(statusQuo = total)
  
  simB <- simB %>%
    bind_rows() %>%
    rename(intervention = total)
  
  # Bind data frames
  compare_df <- left_join(simA, simB) %>%
    mutate(difference = 100 * (intervention - statusQuo) / statusQuo) %>%
    group_by(year, month) %>%
    summarize(mean_difference = mean(difference, na.rm = T),
              sd              = sd(difference),
              upper           = mean_difference + qnorm(0.975) * sd / sqrt(sample_size),
              lower           = mean_difference - qnorm(0.975) * sd / sqrt(sample_size))
  
  # Plot
  compare_plot <- compare_df %>%
    ggplot(aes(x = month, y = mean_difference)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'blue', alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2, color = 'red') +
    scale_y_continuous(labels = comma) +
    labs(x = 'Time (months)',
         y = '% change from status quo',
         title = paste0(metric, ' relative to status quo')) +
    theme_bw()
  
  return(list('compare_plot' = compare_plot, 'compare_df' = compare_df))
}

