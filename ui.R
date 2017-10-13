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
      menuItem("Diagnostics", tabName = 'diagnostics', icon = icon('dashboard')),
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
                    height = 400,
                    p('Note: The model will re-run when these settings are changed, causing model outputs to temporarily fade.'),
                    
                    # Biological inputs
                    sliderInput('depletion', 'Current status relative to unfished population:', 
                                min = 0.1, max = 1, step = 0.1, value = 0.5),
                    
                    # Recruitment month
                    selectInput('recruit_month', 'Recruitment month:', 
                                choices = c(1:12), 
                                selected = 1,
                                multiple = FALSE)),
                
                box(title = 'Biologial Settings', solidHeader = TRUE, status = 'info',
                    width = 4,
                    height = 400,
                    sliderInput('mortality', 'Natural mortality rate:', 
                                min = 0.2, max = 2, step = 0.1, value = 1.76),
                    sliderInput('vbK', 'von Bertalanffy K:', 
                                min = 1, max = 1.5, step = 0.01, value = 1.29),
                    sliderInput('linf', 'L infinity:', 
                                min = 20, max = 24, step = 0.1, value = 22.1)),
                
                box(title = 'Selected Parameters', solidHeader = TRUE, status = 'info',
                    width = 4,
                    height = 400,
                    tableOutput('params'))
                
              )),
      
      # Simulation Model
      tabItem(tabName = "model",
              fluidRow(
                column(width = 3,
                       box(title = 'Simulation Parameters', solidHeader = TRUE, status = 'info',
                           width = NULL,
                           
                           # Recruitment variability
                           selectInput('recruit_vary', 'Set degree of recruitment variability:', 
                                       choices = c('Low' = 0.5, 'Medium' = 1, 'High' = 2),
                                       selected = 1),
                           
                           # Fishing mortality for simulation
                           sliderInput('f_sim', 'Set fishing mortality (F) for simulation:', min = 0, max = 4, step = 0.25, value = 1.75),
                           # Mesh size
                           selectInput('mesh_size', 'Select gillnet mesh size', choices = c('2.54 cm' = 'selA', '3.175 cm' = 'selB'), selected = 2.54),
                           
                           # Closed season months
                           sliderInput('season',
                                       label = 'Set closed season (months):',
                                       min = 1, max = 12, value = c(1,3)),
                           
                           # Length of simulation
                           sliderInput('sim_length', label = 'Length of simulation (months):',
                                       min = 48, max = 120, value = 72))
                ),
                
                column(width = 9,
                       tabBox(title = 'Catch', 
                              width = NULL,
                              tabPanel("Projections",
                                       dygraphOutput('catch', height = '250px')
                              ),
                              tabPanel("Results",
                                       tableOutput('catch_table'))
                       ),
                       tabBox(title = 'Revenue', 
                              width = NULL,
                              tabPanel("Projections",
                                       dygraphOutput('revenue', height = '250px')
                              ),
                              tabPanel("Results",
                                       tableOutput('revenue_table'))
                       ),
                       tabBox(title = 'Biomass', 
                              width = NULL,
                              tabPanel("Projections",
                                       dygraphOutput('sim_biomass', height = '250px')
                              ),
                              tabPanel("Results",
                                       tableOutput('biomass_table'))
                       )
                )
              )),
      
      # Diagnostics
      tabItem(tabName = 'diagnostics',
              fluidRow(
                box(title = 'Fishing Mortality', solidHeader = TRUE, status = 'warning',
                    width = 6),
                box(title = 'Mean Length', solidHeader = TRUE, status = 'warning',
                    width = 6)
              )),
      
      # Methods
      tabItem(tabName = "methods",
              fluidRow(
                box(title = 'Methods', solidHeader = TRUE, status = 'primary',
                    width = 12,
                    includeMarkdown(path = 'sardine_methods_shiny.Rmd'))
              ),
              fluidRow(
                box(title = 'Download', solidHeader = TRUE, status = 'primary',
                    width = 6,
                    # Download reports column
                    downloadButton("methods", "Download Complete Methods"))
              ))
    )
  )
)
