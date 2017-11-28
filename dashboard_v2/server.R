# Sardine App Sever Script -------------------------------------------------
## Project: EDF Sardine Tool
## Script purpose: Server code for Shiny Tool
## Date: June 17th, 2017
## Author: Tyler Clavelle

shinyServer(function(input, output) {
  
  # Functions ------------------------------------------------

  source('../functions/sardine_sim_v3.R')
  source('../functions/sardine_app_functions.R')
  source('../functions/sim_plot.R')
  source('../functions/depletion_v3.R')
  source('../functions/simDepletion.R')
  
  # Catch History Plot ------------------------------------------------------
  
  output$catch_hist <- renderDygraph({
    dygraph(p1) %>%
      dyOptions(fillGraph = TRUE) %>%
      # dySeries(label = "Commercial") %>%
      dySeries(label = "Municipal") %>%
      dySeries(label = "Total") %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      dyRangeSelector(dateWindow = c("2002-04-01", "2017-04-01")) %>%
      dyRibbon(season, palette = c("#ffe6e6","#efefef")) %>%
      dyAxis('y', label = 'Catch (metric tons)') %>%
      dyAxis('x', label = 'Year') 
  })
  
# Baseline Population -----------------------------------------------------

  # Generate length and weight at age and maturity at age vectors given life history parameter settings
  length_weight <- reactive({ 
    # Length at age
    l_at_a <- input$linf * (1-exp(-(1/12)*(input$vbK * (7:max_age - t0))))
    # Weight at age
    w_at_a <- lwA * l_at_a ^ lwB  
    # generate maturity vector
    mat_ramp <- seq(0,1, length.out = 5) # assume proportion mature is linear between ages
    # maturity vector function 
    mat_out <- c(mat_ramp[2:length(mat_ramp)], # maturity ramp
                 rep(1, times = length(l_at_a)-length(mat_ramp) + 1)) # all large fish are mature 
    
    return(list('length_age' = l_at_a, 'weight_age' = w_at_a, 'maturity' = mat_out))
  })
  
  # Initial population settings given intital recruitment settings
  baseline <- reactive({
    # calculate ssb0 from resulting r0
    virgin <- input$r0 * exp(- (1 / 12) * input$mortality * c(0:(max_age - recruit_age)))
    
    # Set recruitment ages. Include max age of analysis
    recruit_ages <- unique(c(seq(from = input$recruit_month, to = max_age, by = 12)))
    
    # Set recruitment months
    recruit_months <- seq(from = as.numeric(input$recruit_month) + 12, to = input$sim_length, by = 12)
    
    # only include numbers for recruitment ages (cohorts)
    virgin[!(c(1:(max_age - recruit_age + 1)) %in% recruit_ages)] <- 0
    
    # set initial pop to virgin
    initial_pop <- virgin
    
    # Get length, weight, and maturity at age vectors
    age_params <- length_weight()
    
    return(list('initial_pop' = initial_pop, 'recruit_months' = recruit_months))
  })
  
  # Recruitment Variability -------------------------------------------------
  
  ### Generate recruitment variability
  recruit_react <- reactive({
    
    # Make rnorm function repeatable
    recruit_rnorm <- repeatable(rlnorm)
    
    # Create an array of recruitment timeseries for running a number of simulations
    rec_var <- array(recruit_rnorm(input$sim_length, log(1), as.numeric(input$recruit_vary)),
                       dim = c(input$sim_number, input$sim_length))

    return(rec_var)
    
  })
  
  # Closed Season Function --------------------------------------------------
  
  # Reactive function to generate the closed season for displaying on the plot
  closed_season <- reactive({
    # generate vector of closed season months 
    closed_period <- seq(from = input$season[1], to = input$season[2], by = 1)
    closed <- sapply(closed_period, function(x) seq(from = x, to = input$sim_length, by = 12)) %>% 
      unlist() %>%
      sort()
    return(closed)
  }) 

# Base Simulation ---------------------------------------------------------
  
  # Run base simulation
  baseModel <- reactive({
    
    # Get length, weight, and maturity at age vectors
    age_params <- length_weight()
    
    # Base simulation
    simBase <- sardine_sim_v3(initial_pop = baseline()[['initial_pop']],
                              recruit_type = 'constant',
                              r0 = input$r0,
                              f = 0,
                              ages = c(7:max_age),
                              sim_length = input$sim_length,
                              length_at_age = age_params[['length_age']],
                              weight_at_age = age_params[['weight_age']],
                              maturity = age_params[['maturity']])
    
    # Update SSB
    sim_m_out <- simBase[['m_out']] * age_params[['weight_age']]
    sim_m_out <- sim_m_out[baseline()[['recruit_months']],]
    
    # Take max of spawning stock biomass from base simulation
    ssb0_wt <- max(rowSums(sim_m_out))
  
    return(list('simBase' = simBase, 'ssb0_wt' = ssb0_wt))
    })

# Run Model ---------------------------------------------------------------

  # Run optimization to find F that matches target depletion
  runModel <- reactive({
    
    # Get length, weight, and maturity at age vectors
    age_params <- length_weight()
    
    # Run baseline model
    baseline <- baseModel()
    
    # Pull out ssb
    ssb0_base <- baseline[['ssb0_wt']]
    
    # Run optimization
    OUT <- optimize(depletion_v3,
                    lower = 0.01,
                    upper = 2,
                    target_depletion = input$depletion,
                    r0 = input$r0, 
                    ssb0_wt = ssb0_base, 
                    initial_pop = baseline()[['initial_pop']], 
                    sim_length = input$sim_length,
                    length_at_age = age_params[['length_age']],
                    weight_at_age = age_params[['weight_age']],
                    maturity = age_params[['maturity']],
                    recruit_months = baseline()[['recruit_months']])
    
    # Run simulation with optimized F to generate initial population for forward simulations
    virgin <- input$r0 * exp(- (1 / 12) * (input$mortality + OUT$minimum) * c(0:(max_age - recruit_age)))
    
    # Set recruitment ages. Include max age of analysis
    recruit_ages <- unique(c(seq(from = input$recruit_month, to = max_age, by = 12)))
    
    # only include numbers for recruitment ages (cohorts)
    virgin[!(c(1:(max_age - recruit_age + 1)) %in% recruit_ages)] <- 0
    
    # set initial pop to resulting
    initial_pop <- virgin
    
    # Run simulation with optimized F to match depletion setting
    simA <- sardine_sim_v3(initial_pop = initial_pop,
                              recruit_type = 'bev_holt',
                              r0 = input$r0,
                              f = OUT$minimum,
                           ssb0 = ssb0_base,
                              ages = c(7:max_age),
                              sim_length = input$sim_length,
                              selectivity = input$mesh_size,
                              recruit_var = recruit_react(),
                              sims = input$sim_number,
                              length_at_age = age_params[['length_age']],
                              weight_at_age = age_params[['weight_age']],
                              maturity = age_params[['maturity']])
    
    # Run simulation with closed season
    simB <- sardine_sim_v3(initial_pop = initial_pop,
                           recruit_type = 'bev_holt',
                           r0 = input$r0,
                           f = OUT$minimum,
                           ssb0 = ssb0_base,
                           ages = c(7:max_age),
                           sim_length = input$sim_length,
                           sims = input$sim_number,
                           selectivity = input$mesh_size,
                           recruit_var = recruit_react(),
                           length_at_age = age_params[['length_age']],
                           weight_at_age = age_params[['weight_age']],
                           maturity = age_params[['maturity']],
                           closed = closed_season())
    
    return(list('optF' = OUT$minimum, 'currentF' = simA, 'simF' = simB, 'ssb_wt' = ssb0_base))
  })
  

# Base Model Plots --------------------------------------------------------
  
  # Base model plots
  baseModelPlots <- reactive({
    
    # Run Model
    simBase <- baseModel()[['simBase']]
  
    # Generate plots
    base_plot_ssb <- sim_plot(sim_df = simBase[['m_out']], sim_length = input$sim_length, metric = 'Spawning Biomass')
    base_plot_catch <- sim_plot(sim_df = simBase[['c_out']], sim_length = input$sim_length, metric = 'Catch')
    
    return(list('ssb' = base_plot_ssb, 'catch' = base_plot_catch))
  })


# Simulation Plots --------------------------------------------------------

  # Simulation model plots
  simModelPlots <- reactive({
    
    # Run model
    sims <- runModel()
    # Get length, weight, and maturity at age vectors
    age_params <- length_weight()
    
    # Extract model runs
    currentF <- sims[['currentF']]
    simF <- sims[['simF']]
    
    # Calculate depletion timeseries
    base_dp <- simDepletion(mature_df = currentF[['m_out']], sim_name = 'Status Quo', 
                            weight_age = age_params[['weight_age']], 
                            sim_length = input$sim_length, ssb0_wt = sims[['ssb_wt']], 
                            recruit_months = baseline()[['recruit_months']])
    
    forward_dp <- simDepletion(mature_df = simF[['m_out']], sim_name = 'Intervention', 
                               weight_age = age_params[['weight_age']], 
                               sim_length = input$sim_length, ssb0_wt = sims[['ssb_wt']], 
                               recruit_months = baseline()[['recruit_months']])
    
    # Combine dataframes
    sim_dp <- bind_rows(base_dp, forward_dp)
    
    # Create depletion plot
    dep_plot_out <- depletion_plot(sim_dp)
    
    # Create catch simulation plot
    catch_plot <- simCompare(simA = currentF[['c_list']], 
                             simB = simF[['c_list']],
                             sim_length = input$sim_length,
                             metric = 'Catch',
                             sample_size = input$sim_number)
    
    # Create biomass simulation plot
    bio_plot <- simCompare(simA = currentF[['m_list']], 
                             simB = simF[['m_list']],
                             sim_length = input$sim_length,
                             metric = 'Biomass',
                           sample_size = input$sim_number)
    
    return(list('depletion_df' = sim_dp, 
                'depletionPlot' = dep_plot_out,
                'catchPlot' = catch_plot[['compare_plot']], 
                'bioPlot' = bio_plot[['compare_plot']]))
  })
  

# Simulation Tables -------------------------------------------------------

  simModelTables <- reactive({
    
    # Run baseline
    base <- baseModel()
    # Run model
    sims <- runModel()
    # Get length, weight, and maturity at age vectors
    age_params <- length_weight()
    
    # Extract model runs
    currentF <- sims[['currentF']]
    simF <- sims[['simF']]
    
    # helper function to process matrices
    goLong <- function(df1, df2, metric) {
      df_out <- df1 %>%
        tbl_df() %>%
        mutate(month     = c(1:nrow(.)),
               year      = ceiling(rep_len(x = c(1:input$sim_length / 12), length.out = input$sim_length + 1)),
               total     = rowSums(.),
               scenario  = 'sq') %>%
        bind_rows(df2 %>%
                    tbl_df() %>%
                    mutate(month     = c(1:nrow(.)),
                           year      = ceiling(rep_len(x = c(1:input$sim_length / 12), length.out = input$sim_length + 1)),
                           total     = rowSums(.),
                           scenario  = 'closed')) %>%
        select(month, year, total, scenario) %>% 
        spread(key = scenario , value = total) %>%
        mutate(metric = metric)
      return(df_out)
    }
    
    # Procees data frames
    catch_df <- goLong(df1 = currentF[['c_out']], df2 = simF[['c_out']], metric = 'Catch')
    bio_df <- goLong(df1 = currentF[['b_out']], df2 = simF[['b_out']], metric = 'Biomass')
    
    # Summarize by year
    catch_table <- catch_df %>%
      group_by(year, metric) %>%
      summarize(sq         = sum(sq, na.rm = T),
                closed     = sum(closed, na.tm = T),
                difference = round(100 * (closed - sq) / sq)) %>%
      select(year, metric, difference)
    
    bio_table <- bio_df %>%
      group_by(year, metric) %>%
      summarize(sq         = sum(sq, na.rm = T),
                closed     = sum(closed, na.tm = T),
                difference = round(100 * (closed - sq) / sq)) %>%
      select(year, metric, difference)
    
    # Summarize over full timeseries for status boxes
    catch_status_box <- catch_df %>%
      summarize(sq         = sum(sq, na.rm = T),
                closed     = sum(closed, na.tm = T),
                difference = round(100 * (closed - sq) / sq)) %>% .$difference
    
    bio_status_box <- bio_df %>%
      summarize(sq         = sum(sq, na.rm = T),
                closed     = sum(closed, na.tm = T),
                difference = round(100 * (closed - sq) / sq)) %>% .$difference
    
    return(list('catch_table' = catch_table, 
                'bio_table' = bio_table,
                'catch_box' = catch_status_box,
                'bio_box'   = bio_status_box))
  })
  

# Render UI Objects -------------------------------------------------------

  # Render plots
  output$base_run_ssb <- renderPlot({ baseModelPlots()[['ssb']] })
  output$base_run_catch <- renderPlot({ baseModelPlots()[['catch']] })
  
  output$sim_depletion <- renderPlot({ simModelPlots()[['depletionPlot']] })
  output$sim_depletion_table <- renderDataTable(simModelPlots()[['depletion_df']])
  
  output$sim_catch_plot <- renderPlot({ simModelPlots()[['catchPlot']] })
  output$sim_bio_plot <- renderPlot({ simModelPlots()[['bioPlot']] })
  
  # Render tables
  output$sim_catch_table <- renderTable({ simModelTables()[['catch_table']]})
  output$sim_bio_table <- renderTable({ simModelTables()[['bio_table']]})
  
  # Render status boxes
  output$catch_status_box <- renderValueBox({
    valueBox(
      simModelTables()[['catch_box']],
      '% Change in catch', icon = icon('ship'), color = 'green')
      })
  
  output$bio_status_box <- renderValueBox({
    valueBox(
      simModelTables()[['bio_box']],
      '% Change in biomass', icon = icon('globe'), color = 'green')
  })

})