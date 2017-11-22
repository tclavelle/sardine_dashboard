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
    OUT <- optimize( 
      depletion_NLL_v2, 
      lower = c(1e6),
      upper = c(1e30),
      depletion = input$depletion,
      catch = catch_history,
      length_at_age = length_at_age,
      weight_at_age = weight_at_age,
      maturity = maturity)
    
    # calculate ssb0 from resulting r0
    virgin <- OUT$minimum[1] * exp(- (1 / 12) * c(0:41))
    
    # Calculate number mature and multiply times weight at age to get SSB0
    ssb0 <- virgin * maturity * weight_at_age
    
    recruit_months <- seq(from = 1, to = 42, by = 12)
    all_months <- c(1:42)
    
    # Only include numbers in recruitment 
    ssb0 <- sum(ssb0[recruit_months])
    
    # Set start population
    virgin[!all_months %in% recruit_months] <- 0
    
    # Simulation with catch history
    fishing_sim <- sardine_optim_v2(start_pop = virgin,
                                    ssb0 = ssb0,
                                    r0 = OUT$minimum[1],
                                    length_at_age = length_at_age,
                                    weight_at_age = weight_at_age,
                                    maturity = maturity,
                                    catch = catch_history,
                                    recruit_type = 'bev_holt')
    return(list(r0 = OUT$minimum[1],
                ssb0 = ssb0,
                history = fishing_sim))
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
  
  # Run Simulation w/ Closed Season -------------------------------------------------
  run_sim <- reactive({
    
    # Run optimization to find population parameters
    pop_params <- run_model()
    
    # Pull out optimized parameters to use in simulations
    r0 <- pop_params[['r0']]
    ssb0 <- pop_params[['ssb0']]
    historic <- pop_params[['history']][['abundance']]
    
    # initial population is last month of historic
    start_pop <- historic[nrow(historic),]
    
    results <- sardine_sim_v2(start_pop = start_pop,
                              ssb0 = ssb0,
                           m = input$mortality,
                           r0 = r0,
                           f = input$f_sim,
                           selectivity = input$mesh_size,
                           length_at_age = length_at_age,
                           weight_at_age = weight_at_age,
                           maturity = maturity,
                           f_mode = 'rate',
                           recruit_type = 'bev_holt',
                           recruit_var = recruit_react(),
                           recruit_month = input$recruit_month,
                           closed = closed_season(),
                           sim_length = input$sim_length)
    
    return(results)
  
    })
  
  # Run Simulation w/ No Closed Season Function -------------------------------------------------
  run_sim_no_closed <- reactive({
    
    # Run optimization to find population parameters
    pop_params <- run_model()
    
    # Pull out optimized parameters to use in simulations
    r0 <- pop_params[['r0']]
    ssb0 <- pop_params[['ssb0']]
    historic <- pop_params[['history']][['abundance']]
    
    # initial population is last month of historic
    start_pop <- historic[nrow(historic),]
    
    results <- sardine_sim_v2(start_pop = start_pop,
                              ssb0 = ssb0,
                              m = input$mortality,
                              r0 = r0,
                              f = input$f_sim,
                              selectivity = input$mesh_size,
                              length_at_age = length_at_age,
                              weight_at_age = weight_at_age,
                              maturity = maturity,
                              f_mode = 'rate',
                              recruit_type = 'bev_holt',
                              recruit_var = recruit_react(),
                              recruit_month = input$recruit_month,
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
    opt_params <- data_frame(Parameter = c('Virgin SSB (MMT)',
                                           'Initial recruitment (millions)'),
                             value = c(pop_params[['ssb0']] / 1e6,
                                       pop_params[['r0']] / 1e6))
    return(opt_params)
  
    })
  
  # Summary results tables
  catch_results_summary <- reactive({

    # Catch
    sim_out_catch <- run_sim()[['c_out']] %>%
      tbl_df() %>%
      mutate(month = rep_len(1:12, length.out = input$sim_length + 1),
             year = ceiling(rep_len(x = c(1:input$sim_length / 12), length.out = input$sim_length + 1))) %>%
      gather(key = 'age', value = 'number', 1:c(ncol(.)-2)) %>%
      group_by(year) %>%
      summarize(Value = sum(number, na.rm = T) / 1e5) %>% 
      mutate(Metric = 'Catch (100,000s MT)',
             Policy = 'Closed Season') %>%
      spread(Policy, Value)
    
    sim_out_no_closed_catch <- run_sim_no_closed()[['c_out']] %>%
      tbl_df() %>%
      mutate(month = rep_len(1:12, length.out = input$sim_length + 1),
             year = ceiling(rep_len(x = c(1:input$sim_length / 12), length.out = input$sim_length + 1))) %>%
      gather(key = 'age', value = 'number', 1:c(ncol(.)-1)) %>%
      group_by(year) %>%
      summarize(Value = sum(number, na.rm = T) / 1e5) %>% 
      mutate(Metric = 'Catch (100,000s MT)',
             Policy = 'No Closed Season') %>%
      spread(Policy, Value)
    
    # Join Catch Tables
    sim_out_catch <- left_join(sim_out_catch, sim_out_no_closed_catch) %>%
      mutate(Difference = `Closed Season` - `No Closed Season`)
    
    return(sim_out_catch)
  
    })
  
  biomass_results_summary <- reactive({
    # Biomass
    sim_out_biomass <- run_sim()[['b_out']] %>%
      tbl_df() %>%
      mutate(month = rep_len(1:12, length.out = input$sim_length + 1),
             year = ceiling(rep_len(x = c(1:input$sim_length / 12), length.out = input$sim_length + 1))) %>%
      gather(key = 'age', value = 'number', 1:c(ncol(.)-2)) %>%
      filter(month == 12) %>%
      group_by(year) %>%
      summarize(Value = sum(number, na.rm = T)) %>% 
      mutate(Metric = 'Biomass (MT)',
             Policy = 'Closed Season') %>%
      spread(Policy, Value)
    
    sim_out_no_closed_biomass <- run_sim_no_closed()[['b_out']] %>%
      tbl_df() %>%
      mutate(month = rep_len(1:12, length.out = input$sim_length + 1),
             year = ceiling(rep_len(x = c(1:input$sim_length / 12), length.out = input$sim_length + 1))) %>%
      gather(key = 'age', value = 'number', 1:c(ncol(.)-2)) %>%
      filter(month == 12) %>%
      group_by(year) %>%
      summarize(Value = sum(number, na.rm = T)) %>% 
      mutate(Metric = 'Biomass (MT)',
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

output$biomass_table <- renderTable({
    # Run summary function
    biomass_results_summary()
  })
  
  # Result Plots -------------------------------------------------
  
  ## Equilibrium Biomass
  # output$biomass <- renderDygraph({
  #   
  #   base <- run_model()[['base']]
  #   
  #   model_out_summary <- base[[2]] %>%
  #     tbl_df() %>%
  #     group_by(month) %>%
  #     summarize(biomass   = sum(biomass, na.rm = T) / 1000)
  #   
  #   dygraph(model_out_summary) %>%
  #     dyAxis('y', label = 'Biomass (1000s MT)') %>%
  #     dyAxis('x', label = 'Month') %>%
  #     dyOptions(colors = RColorBrewer::brewer.pal(2, "Dark2"))
  #   
  # })
  # 
  ## Historical Biomass
  # # output$hist_biomass <- renderDygraph({
  #   
  #   history <- run_model()[['history']]
  #   
  #   model_out_summary <- history[['biomass']] %>%
  #     tbl_df() %>%
  #     mutate(month = c(1:157),
  #            year  = ceiling(rep_len(x = c(1:64 / 12), length.out = 157))) %>%
  #     select(month, year, everything()) %>%
  #     gather(key = 'age_class', value = 'biomass', 3:ncol(.)) %>%
  #     group_by(month) %>%
  #     summarize(Biomass = sum(biomass, na.rm = T) / 1000)
  #   
  #   dygraph(model_out_summary) %>%
  #     dyAxis('y', label = 'Biomass (1000s MT)') %>%
  #     dyAxis('x', label = 'Month') %>%
  #     dyOptions(colors = RColorBrewer::brewer.pal(2, "Dark2"))
  #   
  # })
  # 
  ## Projected Catch   
  output$catch <- renderDygraph({

    sim_out_summary <- run_sim()[['c_out']] %>%
      tbl_df() %>%
      mutate(month = c(1:nrow(.))) %>%
      gather(key = 'age', value = 'number', 1:c(ncol(.)-1)) %>%
      group_by(month) %>%
      summarize(`Catch` = sum(number, na.rm = T) / 1000)

    sim_out_summary_no_closed <- run_sim_no_closed()[['c_out']] %>%
      tbl_df() %>%
      mutate(month = c(1:nrow(.))) %>%
      gather(key = 'age', value = 'number', 1:c(ncol(.)-1)) %>%
      group_by(month) %>%
      summarize(`Catch (no closed season)` = sum(number, na.rm = T) / 1000)

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

  ## Projected Biomass  
  output$sim_biomass <- renderDygraph({
    
    # Run simulation
    sim_out_summary <- run_sim()[['b_out']] %>%
      tbl_df() %>%
      mutate(month = c(1:nrow(.))) %>%
      gather(key = 'age', value = 'number', 1:c(ncol(.)-1)) %>%
      group_by(month) %>%
      summarize(Biomass = sum(number, na.rm = T) / 1000)
    
    sim_out_summary_no_closed <- run_sim_no_closed()[['b_out']] %>%
      tbl_df() %>%
      mutate(month = c(1:nrow(.))) %>%
      gather(key = 'age', value = 'number', 1:c(ncol(.)-1)) %>%
      group_by(month) %>%
      summarize(`Biomass (no closed season)` = sum(number, na.rm = T) / 1000)
    
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