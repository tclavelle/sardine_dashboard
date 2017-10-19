# Sardine App Sever Script -------------------------------------------------
## Project: EDF Sardine Tool
## Script purpose: Server code for Shiny Tool
## Date: June 17th, 2017
## Author: Tyler Clavelle

shinyServer(function(input, output) {
  
  # Basic Fishery Demography ------------------------------------------------
  
  vbK <- 1.29 # von Bert K (1.29)
  Linf <- 22.1 # L infinity (22.1)
  lwA <- 0.0078 # alpha from length-weight relationship
  lwB <- 3.165 # beta from length weight relationship
  t0 <- -lwA/lwB
  M <- 1.786 # natural mortality (1.786)
  set.seed(100)
  
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
  
  # Run Model Function ------------------------------------------------------
  run_model <- reactive({
    
    # Run optimization to find initial recruitment 
    OUT <- optim(c(1e4, 1), 
                 depletion_NLL, 
                 method = 'L-BFGS-B',
                 lower = c(1e2, 0),
                 upper = c(1e10, 2),
                 depletion = input$depletion,
                 catch = catch_history,
                 soi = exp(soi_monthly$SOI))
    
    # Base simulation
    base_sim <- sardine_sim(r0 = OUT$par[1], 
                            m = input$mortality,
                            linf = input$linf,
                            vbk = input$vbK,
                            f=0, 
                            recruit_month = input$recruit_month, 
                            recruit_type = 'constant')
    
    # Set B0 to sum of spawning stock biomass in numbers ()
    B0 <- sum(base_sim[[6]], na.rm = T)
    ssb0 <- B0
    
    # Optimized value of R0
    R0 <-  OUT$par[1]
    
    # Assume steepness of 0.8
    h  <- 0.8
    
    # Beverton-Holt parameters
    alpha_bh		<-	(B0 / R0) * (1-h)/(4*h)
    beta_bh		<-	(5*h-1)/(4*h*R0)
    
    # Adjust SOI by SOI parameter value
    soi_for_opt <- exp(soi_monthly$SOI) * OUT$par[2]
    
    # Simulation with catch history
    fishing_sim <- sardine_optim(ssb0 = ssb0,
                                 r0 = R0,
                                 catch = catch_history, 
                                 soi = soi_for_opt, 
                                 recruit_month = 1, 
                                 alpha_bh = alpha_bh, 
                                 beta_bh = beta_bh, 
                                 recruit_type = 'bev_holt')
    # browser()
    return(list(r0 = OUT$par[1],
                alpha_bh = alpha_bh,
                beta_bh = beta_bh,
                base = base_sim,
                ssb0 = ssb0,
                history = fishing_sim,
                soi_param = OUT$par[2]))
  })
  
  # Recruitment Variability -------------------------------------------------
  
  ### Generate recruitment variability
  recruit_react <- reactive({
    
    # Make rnorm function repeatable
    recruit_rnorm <- repeatable(rnorm)
    # Set recruitment variability (may need to make this its own isolated function for Shiny)
    rec_var	<-	lapply(input$sim_length, FUN = recruit_rnorm, sd = as.numeric(input$recruit_vary)) %>%
      unlist()
    rec_var <- exp(rec_var)
    return(rec_var)
    
    })
  
  ### SOI parameter for projection
  soi_react <- reactive({
    # Calculate mean and SD of SOI timeseries
    soi_mean <- mean(soi$SOI, na.rm = T)
    soi_sd <- sd(soi$SOI, na.rm = T)
    # Make sample repeatable
    # soi_sample <- repeatable(sample)
    soi_rnorm <- repeatable(rnorm)
    # Sample, with replacement, soi values for duration of projection
    # soi_project <- soi_sample(x = soi$SOI, size = input$sim_length, replace = TRUE)
    soi_project	<-	lapply(input$sim_length, FUN = soi_rnorm, mean = soi_mean, sd = soi_sd) %>%
      unlist()
    # Exponentiate values to convert negatives to decimals
    soi_project <- exp(soi_project)
    return(soi_project)
    
    })
  
  # Run Simulation w/ Closed Season -------------------------------------------------
  run_sim <- reactive({
    
    # Run optimization to find population parameters
    pop_params <- run_model()
    
    # Pull out optimized parameters to use in simulations
    r0 <- pop_params[['r0']]
    alpha <- pop_params[['alpha_bh']]
    beta <- pop_params[['beta_bh']]
    soi_param <- pop_params[['soi_param']] # soi variability parameter
    
    # Adjust SOI by the optimized parameter
    soi_project <- soi_react() * soi_param
    
    results <- sardine_sim(m = input$mortality,
                           r0 = r0,
                           vbk = input$vbK,
                           linf = input$linf,
                           alpha_bh = alpha,
                           beta_bh = beta,
                           f = input$f_sim,
                           selectivity = input$mesh_size,
                           f_mode = 'rate',
                           recruit_type = 'bev_holt',
                           recruit_var = recruit_react(),
                           soi_project = soi_project,
                           recruit_month = input$recruit_month,
                           closed = closed_season(),
                           sim_length = input$sim_length)
    
    return(results)
  
    })
  
  # Run Simulation w/ No Closed Season Function -------------------------------------------------
  run_sim_no_closed <- reactive({
    
    # Run optimization to find population parameters
    pop_params <- run_model()
    
    r0 <- pop_params[['r0']]
    alpha <- pop_params[['alpha_bh']]
    beta <- pop_params[['beta_bh']]
    soi_param <- pop_params[['soi_param']] # soi variability parameter

    # Adjust SOI by the optimized parameter
    soi_project <- soi_react() * soi_param
    
    results <- sardine_sim(m = input$mortality,
                           r0 = r0,
                           vbk = input$vbK,
                           linf = input$linf,
                           alpha_bh = alpha,
                           beta_bh = beta,
                           f = input$f_sim,
                           selectivity = input$mesh_size,
                           f_mode = 'rate',
                           recruit_type = 'bev_holt',
                           recruit_var = recruit_react(),
                           soi_project = soi_project,
                           closed = NA,
                           sim_length = input$sim_length)
    
    return(results)
  
    })
  
  # Result Tables -------------------------------------------------
  
  # Biology Parameter Table
  output$params <- renderTable({
    
    params <- data_frame(Parameter = c('Natural mortality',
                                       'von Bertalanffy K', 
                                       'L infinity',
                                       'Current depletion',
                                       'Recruitment month'),
                         Value = c(input$mortality,
                                   input$vbK,
                                   input$linf,
                                   input$depletion,
                                   input$recruit_month))
  
    })
  
  # Optimized parameter table
  output$opt_param_table <- renderTable({
    
    # Run model
    pop_params <- run_model()
    
    # Make table
    opt_params <- data_frame(Parameter = c('Virgin SSB',
                                           'Initial recruitment',
                                           'Alpha', 
                                           'Beta'),
                             value = c(pop_params[['ssb0']],
                                       pop_params[['r0']],
                                       pop_params[['alpha_bh']],
                                       pop_params[['beta_bh']]))
    return(opt_params)
  
    })
  
  # Summary results tables
  catch_results_summary <- reactive({
    
    # Catch
    sim_out_catch <- run_sim()[[3]] %>%
      tbl_df() %>%
      group_by(year) %>%
      summarize(Value = sum(catch, na.rm = T) / 1000) %>% 
      mutate(Metric = 'Catch (1000s MT)',
             Policy = 'Closed Season') %>%
      spread(Policy, Value)
    
    sim_out_no_closed_catch <- run_sim_no_closed()[[3]] %>%
      tbl_df() %>%
      group_by(year) %>%
      summarize(Value = sum(catch, na.rm = T) / 1000) %>%
      mutate(Metric = 'Catch (1000s MT)',
             Policy = 'No Closed Season') %>%
      spread(Policy, Value)
    
    # Join Catch Tables
    sim_out_catch <- left_join(sim_out_catch, sim_out_no_closed_catch) %>%
      mutate(Difference = `Closed Season` - `No Closed Season`)
    
    return(sim_out_catch)
  
    })
  
  
  revenue_results_summary <- reactive({
    
    # Revenue
    sim_out_revenue <- run_sim()[[3]] %>%
      tbl_df() %>%
      group_by(year) %>%
      summarize(Value = sum(revenue, na.rm = T) / 1e6) %>%
      mutate(Metric = 'Revenue (millions)',
             Policy = 'Closed Season') %>%
      spread(Policy, Value)
    
    sim_out_no_closed_revenue <- run_sim_no_closed()[[3]] %>%
      tbl_df() %>%
      group_by(year) %>%
      summarize(Value = sum(revenue, na.rm = T) / 1e6) %>%
      mutate(Metric = 'Revenue (millions)',
             Policy = 'No Closed Season') %>%
      spread(Policy, Value)
    
    # Join Revenue Result tables
    sim_out_revenue <- left_join(sim_out_revenue, sim_out_no_closed_revenue) %>%
      mutate(Difference = `Closed Season` - `No Closed Season`)
    
    return(sim_out_revenue)
  })
  
  biomass_results_summary <- reactive({
    # Biomass
    sim_out_biomass <- run_sim()[[2]] %>%
      tbl_df() %>%
      group_by(year) %>%
      filter(month == max(month, na.rm = T)) %>%
      summarize(Value = sum(biomass, na.rm = T) / 1000) %>%
      mutate(Metric = 'Biomass (1000s MT)',
             Policy = 'Closed Season') %>%
      spread(Policy, Value)
    
    sim_out_no_closed_biomass <- run_sim_no_closed()[[2]] %>%
      tbl_df() %>%
      group_by(year) %>%
      filter(month == max(month, na.rm = T)) %>%
      summarize(Value = sum(biomass, na.rm = T) / 1000) %>%
      mutate(Metric = 'Biomass (1000s MT)',
             Policy = 'No Closed Season') %>%
      spread(Policy, Value)
    
    # Join Biomass results
    sim_out_biomass <- left_join(sim_out_biomass, sim_out_no_closed_biomass) %>%
      mutate(Difference = `Closed Season` - `No Closed Season`)
    
    return(sim_out_biomass)
  })
  
  output$catch_table <- renderTable({
    # Run summary function
    catch_results_summary()
  })
  output$revenue_table <- renderTable({
    # Run summary function
    revenue_results_summary()
  })
  output$biomass_table <- renderTable({
    # Run summary function
    biomass_results_summary()
  })
  
  # Result Plots -------------------------------------------------
  
  ## Equilibrium Biomass
  output$biomass <- renderDygraph({
    
    base <- run_model()[['base']]
    
    model_out_summary <- base[[2]] %>%
      tbl_df() %>%
      group_by(month) %>%
      summarize(biomass   = sum(biomass, na.rm = T) / 1000)
    
    dygraph(model_out_summary) %>%
      dyAxis('y', label = 'Biomass (1000s MT)') %>%
      dyAxis('x', label = 'Month') %>%
      dyOptions(colors = RColorBrewer::brewer.pal(2, "Dark2"))
    
  })
  
  ## Historical Biomass
  output$hist_biomass <- renderDygraph({
    
    history <- run_model()[['history']]
    
    model_out_summary <- history[['b_out']] %>%
      tbl_df() %>%
      mutate(month = c(1:(64)),
             year  = ceiling(rep_len(x = c(1:64 /12), length.out = 64))) %>%
      select(month, year, everything()) %>%
      gather(key = 'age_class', value = 'biomass', 3:ncol(.)) %>%
      group_by(month) %>%
      summarize(Biomass = sum(biomass, na.rm = T) / 1000)
    
    dygraph(model_out_summary) %>%
      dyAxis('y', label = 'Biomass (1000s MT)') %>%
      dyAxis('x', label = 'Month') %>%
      dyOptions(colors = RColorBrewer::brewer.pal(2, "Dark2"))
    
  })
  
  ## Projected Catch   
  output$catch <- renderDygraph({
    
    sim_out_summary <- run_sim()[[3]] %>%
      tbl_df() %>%
      group_by(month) %>%
      summarize(`Catch` = sum(catch, na.rm = T) / 1000)
    
    sim_out_summary_no_closed <- run_sim_no_closed()[[3]] %>%
      tbl_df() %>%
      group_by(month) %>%
      summarize(`Catch (no closed season)` = sum(catch, na.rm = T) / 1000)
    
    sim_out_summary <- sim_out_summary %>%
      left_join(sim_out_summary_no_closed)
    
    # Build closed season matrix to add to plot
    closed <- closed_season() %>%
      seqToIntervals()
    
    # Base plot
    dy_catch <- dygraph(sim_out_summary,
                        group = 'projections') %>%
      dyAxis('y', label = 'Catch (1000s MT)') %>%
      dyAxis('x', label = 'Month') %>%
      dyLegend(width = 400) %>%
      dyOptions(colors = RColorBrewer::brewer.pal(2, "Dark2"))
    
    # Add shading for each closed season
    for(a in 1:nrow(closed)) {
      dy_catch <- dy_catch %>%
        dyShading(from = closed[a,1], to = closed[a,2], color = '#eba9a9')
    }
    # print plot
    dy_catch
  })  
  
  ## Projected difference in catch   
  output$catch_difference <- renderPlot({
    
    sim_out_summary <- run_sim()[[3]] %>%
      tbl_df() %>%
      group_by(month) %>%
      summarize(`Catch` = sum(catch, na.rm = T) / 1000)
    
    sim_out_summary_no_closed <- run_sim_no_closed()[[3]] %>%
      tbl_df() %>%
      group_by(month) %>%
      summarize(`Catch (no closed season)` = sum(catch, na.rm = T) / 1000)
    
    sim_out_summary <- sim_out_summary %>%
      left_join(sim_out_summary_no_closed) %>%
      mutate(Difference = Catch - `Catch (no closed season)`,
             Sign = ifelse(Difference >= 0, 'pos','neg')) %>%
      select(month, Difference, Sign)
    
    ggplot(sim_out_summary, aes(x = month, y = Difference, fill = Sign)) +
      geom_bar(stat = 'identity') +
      scale_fill_manual(values = c('pos' = '#228B22', 'neg' = '#FF0000')) +
      guides(fill = F) +
      theme_bw()
  }) 
  
  ## Projected Revenue   
  output$revenue <- renderDygraph({
    
    sim_out_summary <- run_sim()[[3]] %>%
      tbl_df() %>%
      group_by(month) %>%
      summarize(`Revenue` = sum(revenue, na.rm = T) / 1e6)
    
    sim_out_summary_no_closed <- run_sim_no_closed()[[3]] %>%
      tbl_df() %>%
      group_by(month) %>%
      summarize(`Revenue (no closed season)` = sum(revenue, na.rm = T) / 1e6)
    
    sim_out_summary <- sim_out_summary %>%
      left_join(sim_out_summary_no_closed)
    
    # Build closed season matrix to add to plot
    closed <- closed_season() %>%
      seqToIntervals()
    
    # Base plot
    dy_revenue <- dygraph(sim_out_summary,
                          group = 'projections') %>%
      dyAxis('y', label = 'Revenue (Millions, $USD)') %>%
      dyAxis('x', label = 'Month') %>%
      dyLegend(width = 400) %>%
      dyOptions(colors = RColorBrewer::brewer.pal(2, "Dark2"))
    
    # Add shading for each closed season
    for(a in 1:nrow(closed)) {
      dy_revenue <- dy_revenue %>%
        dyShading(from = closed[a,1], to = closed[a,2], color = '#eba9a9')
    }
    # print plot
    dy_revenue
  })  
  
  ## Projected Biomass  
  output$sim_biomass <- renderDygraph({
    
    # Run simulation
    sim_out_summary <- run_sim()[[2]] %>%
      tbl_df() %>%
      group_by(month) %>%
      summarize(Biomass = sum(biomass, na.rm = T) / 1000)
    
    sim_out_summary_no_closed <- run_sim_no_closed()[[2]] %>%
      tbl_df() %>%
      group_by(month) %>%
      summarize(`Biomass (no closed season)` = sum(biomass, na.rm = T) / 1000)
    
    sim_out_summary <- sim_out_summary %>%
      left_join(sim_out_summary_no_closed)
    
    # Build closed season matrix to add to plot
    closed <- closed_season() %>%
      seqToIntervals()
    
    # Base plot
    dy_bio <- dygraph(sim_out_summary,
                      group = 'projections') %>%
      dyAxis('y', label = 'Biomass (1000s MT)') %>%
      dyAxis('x', label = 'Month') %>%
      dyLegend(width = 400) %>%
      dyOptions(colors = RColorBrewer::brewer.pal(2, "Dark2"))
    
    # Add shading for each closed season
    for(a in 1:nrow(closed)) {
      dy_bio <- dy_bio %>%
        dyShading(from = closed[a,1], to = closed[a,2], color = '#eba9a9')
    }
    # print plot
    dy_bio
  })
  
  # Methods & Downloads -------------------------------------------------
  
  # Download Methods
  output$methods <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "sardine_methods.html",
    content = function(file) {
      rmarkdown::render('sardine_app_methods.Rmd', output_file = file)
    }
  )
  
})