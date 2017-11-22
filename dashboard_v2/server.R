# Sardine App Sever Script -------------------------------------------------
## Project: EDF Sardine Tool
## Script purpose: Server code for Shiny Tool
## Date: June 17th, 2017
## Author: Tyler Clavelle

shinyServer(function(input, output) {
  
  # Basic Fishery Demography ------------------------------------------------

  source('../functions/sardine_sim_v3.R')
  source('../functions/sardine_app_functions.R')
  source('../functions/sim_plot.R')
  source('../functions/depletion_v3.R')
  source('../functions/simDepletion.R')
 
  lwA <- 0.0078 # alpha from length-weight relationship
  lwB <- 3.165 # beta from length weight relationship
  t0 <- -lwA/lwB
  max_age <- 31
  recruit_age <- 7
  ages <- c(7:max_age)

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
    recruit_months <- seq(from = recruit_month + 12, to = input$sim_length, by = 12)
    
    # only include numbers for recruitment ages (cohorts)
    virgin[!(c(1:(max_age - recruit_age + 1)) %in% recruit_ages)] <- 0
    
    # set initial pop to virgin
    initial_pop <- virgin
    
    # Get length, weight, and maturity at age vectors
    age_params <- length_weight()
    
    return(list('initial_pop' = initial_pop, 'recruit_months' = recruit_months))
  })

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
    ssb0 <- sim_m_out[which(rowSums(sim_m_out)==ssb0_wt),] / age_params[['weight_age']]
  
    return(list('simBase' = simBase, 'ssb0_wt' = ssb0_wt))
    })
  
  # Run optimization to find F that matches target depletion
  runModel <- reactive({
    
    # Get length, weight, and maturity at age vectors
    age_params <- length_weight()
    
    # Run baseline model
    baseline <- baseModel()
    
    # Pull out ssb
    ssb0_wt <- baseline[['ssb0_wt']]
    
    # Run optimization
    OUT <- optimize(depletion_v3,
                    lower = 0.01,
                    upper = 2,
                    target_depletion = input$depletion,
                    r0 = input$r0, 
                    ssb0_wt = ssb0_wt, 
                    initial_pop = baseline()[['initial_pop']], 
                    sim_length = input$sim_length,
                    length_at_age = age_params[['length_age']],
                    weight_at_age = age_params[['weight_age']])
    
    # Run simulation with optimized F to match depletion setting
    simA <- sardine_sim_v3(initial_pop = baseline()[['initial_pop']],
                              recruit_type = 'bev_holt',
                              r0 = input$r0,
                              f = OUT$minimum,
                              ages = c(7:max_age),
                              sim_length = input$sim_length,
                              length_at_age = age_params[['length_age']],
                              weight_at_age = age_params[['weight_age']],
                              maturity = age_params[['maturity']])
    
    # Run simulation with F setting
    simB <- sardine_sim_v3(initial_pop = baseline()[['initial_pop']],
                           recruit_type = 'bev_holt',
                           r0 = input$r0,
                           f = input$f_sim,
                           ages = c(7:max_age),
                           sim_length = input$sim_length,
                           length_at_age = age_params[['length_age']],
                           weight_at_age = age_params[['weight_age']],
                           maturity = age_params[['maturity']])
    
    return(list('optF' = OUT$minimum, 'currentF' = simA, 'simF' = simB, 'ssb_wt' = ssb0_wt))
  })
  
  # Base model plots
  baseModelPlots <- reactive({
    
    # Run Model
    simBase <- baseModel()[['simBase']]
  
    # Generate plots
    base_plot_ssb <- sim_plot(sim_df = simBase[['m_out']], sim_length = input$sim_length, metric = 'Spawning Biomass')
    base_plot_catch <- sim_plot(sim_df = simBase[['c_out']], sim_length = input$sim_length, metric = 'Catch')
    
    return(list('ssb' = base_plot_ssb, 'catch' = base_plot_catch))
  })
  
  # Simulation model plots
  simModelPlots <- reactive({
    
    # Run baseline
    base <- baseModel()
    
    # Run model
    sims <- runModel()
    
    # Get length, weight, and maturity at age vectors
    age_params <- length_weight()
    
    # Extract model runs
    currentF <- sims[['currentF']]
    simF <- sims[['simF']]
    
    # Calculate depletion timeseries
    base_dp <- simDepletion(mature_df = currentF[['m_out']], sim_name = 'current F', 
                            weight_age = age_params[['weight_age']], 
                            sim_length = input$sim_length, ssb0 = sims[['ssb_wt']], 
                            recruit_months = baseline()[['recruit_months']])
    
    forward_dp <- simDepletion(mature_df = simF[['m_out']], sim_name = 'selected F', 
                               weight_age = age_params[['weight_age']], 
                               sim_length = input$sim_length, ssb0 = sims[['ssb_wt']], 
                               recruit_months = baseline()[['recruit_months']])
    
    # Combine dataframes
    sim_dp <- bind_rows(base_dp, forward_dp)
    
    # Create depletion plot
    dep_plot_out <- depletion_plot(sim_dp)
    
    # Create catch simulation plot
    catch_plot <- simCompare(simA = currentF[['c_out']], 
                             simB = simF[['c_out']],
                             sim_length = input$sim_length,
                             metric = 'Catch')
    
    # Create biomass simulation plot
    bio_plot <- simCompare(simA = currentF[['b_out']], 
                             simB = simF[['b_out']],
                             sim_length = input$sim_length,
                             metric = 'Biomass')
    
    return(list('depletion_df' = sim_dp, 'depletionPlot' = dep_plot_out,
                'catchPlot' = catch_plot, 'bioPlot' = bio_plot))
  })
  
  # Render plots
  output$base_run_ssb <- renderPlot({ baseModelPlots()[['ssb']] })
  output$base_run_catch <- renderPlot({ baseModelPlots()[['catch']] })
  
  output$sim_depletion <- renderPlot({ simModelPlots()[['depletionPlot']] })
  output$sim_depletion_table <- renderDataTable(simModelPlots()[['depletion_df']])
  
  output$sim_catch_plot <- renderPlot({ simModelPlots()[['catchPlot']] })
  output$sim_bio_plot <- renderPlot({ simModelPlots()[['bioPlot']] })
  
  # Depletion value !! just for testing
  output$optF <- renderText(runModel()[['optF']])


})