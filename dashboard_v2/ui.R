## ui.R ##
library(shinydashboard)

dashboardPage(
  dashboardHeader(title = 'Fishery Simulator'),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("home")),
      menuItem("The Fishery", tabName = "fishery", icon = icon("ship")),
      menuItem('Parameters', tabName = 'parameters', icon = icon('cog')),
      menuItem("Fishery Model", tabName = "model", icon = icon("area-chart")),
      menuItem("Methods", tabName = "methods", icon = icon("list"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Overview
      tabItem(tabName = "overview",
              fluidRow(
                box(title = 'Overview', solidHeader = TRUE, status = 'primary',
                    width = 12,
                    includeMarkdown(path = 'sardine_overview.Rmd'))
              )
      ), 
      
      # The Fishery
      tabItem(tabName = "fishery",
              fluidRow(
                box(title = 'Landings', solidHeader = TRUE, status = 'primary',
                    width = 12,
                    dygraphOutput('catch_hist'))),
              fluidRow(
                box(title = 'Background', solidHeader = TRUE, status = 'primary',
                    width = 12,
                    includeMarkdown(path = 'sardine_fishery.Rmd'))
              ) 
      ),
      
      # Population Parameters
      tabItem(tabName = 'parameters',
              fluidRow(
                box(title = 'Population Settings', solidHeader = TRUE, status = 'info',
                    width = 4,
                    
                    # Biological inputs
                    sliderInput('depletion', 'Current status relative to unfished population:', 
                                min = 0.1, max = 1, step = 0.1, value = 0.5),
                    
                    sliderInput('mortality', 'Natural mortality rate:', 
                                min = 0.2, max = 2, step = 0.1, value = 1.76),
                    sliderInput('vbK', 'von Bertalanffy K:', 
                                min = 1, max = 1.5, step = 0.01, value = 1.29),
                    sliderInput('linf', 'L infinity:', 
                                min = 20, max = 24, step = 0.1, value = 22.1)),
                
                box(title = 'Simulation Settings', solidHeader = TRUE, status = 'info',
                    width = 4,
                    
                    # Recruitment 
                    numericInput('r0', 'Initial recruitment (individuals)',
                                 value = 1e6),
                    
                    selectInput('recruit_month', 'Recruitment month:', 
                                choices = c(1:12), 
                                selected = 1,
                                multiple = FALSE),
                    
                    # Recruitment variability
                    selectInput('recruit_vary', 'Set degree of recruitment variability:', 
                                choices = c('Low' = 0.1, 'Medium' = 0.25, 'High' = 0.5),
                                selected = 1)#,
                    
                    # Fishing mortality for simulation
                    # sliderInput('f_sim', 'Set fishing mortality (F) for simulation:', min = 0, max = 4, step = 0.25, value = 1.75)
                ),
                
                box(title = 'Model Settings', solidHeader = TRUE, status = 'info',
                    width = 4,
                    
                    # Length of simulation
                    sliderInput('sim_length', label = 'Length of simulation (months):',
                                min = 48, max = 120, value = 72),
                    
                    # Number of simulations
                    sliderInput('sim_number', label = 'Number of simulations to run:',
                                min = 1, max = 100, value = 25)
              )
              )       
      ),
      
      # Simulation Model
      tabItem(tabName = "model",
              fluidRow(
                column(width = 3,
                       box(solidHeader = TRUE, status = 'info',
                           width = NA,
                           # Mesh size
                           selectInput('mesh_size', 'Select gillnet mesh size', choices = c('2.54 cm' = 'selA', '3.175 cm' = 'selB'), selected = 2.54)
                       ) 
                ),
                
                column(width = 3,
                       box(solidHeader = TRUE, status = 'info',
                           width = NA,
                           # Closed season months
                           sliderInput('season',
                                       label = 'Set closed season (months):',
                                       min = 1, max = 12, value = c(1,3))
                       )
                ), 
                
                column(width = 3,
                       valueBoxOutput('catch_status_box', width = NA)
                ),
                
                column(width = 3, 
                       valueBoxOutput('bio_status_box', width = NA)
                )
              ),
              
              fluidRow(
                tabBox(title = 'Catch', 
                       width = NULL,
                       # height = 250,
                       tabPanel("Projections",
                                plotOutput('sim_catch_plot', height = '250px')),
                       tabPanel("Results",
                       tableOutput('sim_catch_table'))
                ),
                tabBox(title = 'Biomass', 
                       width = NULL,
                       # height = 250, 
                       tabPanel("Projections",
                                plotOutput('sim_bio_plot', height = '250px')),
                       tabPanel("Results",
                                tableOutput('sim_bio_table'))
                       
                ),
                tabBox(title = 'Depletion', 
                       width = NULL,
                       # height = 250,
                       tabPanel("Projections",
                                plotOutput('sim_depletion', height = '250px')), 
                       tabPanel("Results",
                                tableOutput('depletion_table'))
                )
              )
      ),
      
      # Methods
      tabItem(tabName = "methods",
              fluidRow(
                box(title = 'Methods Overview', solidHeader = TRUE, status = 'primary',
                    width = 12,
                    includeMarkdown(path = 'sardine_methods_shiny.Rmd'))
              ),
              fluidRow(
                box(title = 'Biological and Fishery Parameters', solidHeader = TRUE, status = 'primary',
                    width = 12, collapsible = TRUE, collapsed = TRUE,
                    withMathJax(includeMarkdown(path = 'sardine_methods_params.Rmd')))
              ),
              fluidRow(
                box(title = 'Biological Model', solidHeader = TRUE, status = 'primary',
                    width = 12, collapsible = TRUE, collapsed = TRUE,
                    withMathJax(includeMarkdown(path = 'sardine_methods_biomodel.Rmd')))
              ),
              fluidRow(
                box(title = 'Current Depletion and Fishing Mortality', solidHeader = TRUE, status = 'primary',
                    width = 12, collapsible = TRUE, collapsed = TRUE,
                    withMathJax(includeMarkdown(path = 'sardine_methods_depletion.Rmd')))
              ),
              fluidRow(
                box(title = 'Policy Simulation', solidHeader = TRUE, status = 'primary',
                    width = 12, collapsible = TRUE, collapsed = TRUE,
                    withMathJax(includeMarkdown(path = 'sardine_methods_policies.Rmd')))
              )
              )
    ) 
  )
)
