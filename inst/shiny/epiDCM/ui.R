shinyUI(pageWithSidebar(
  
  # Header
  headerPanel('EpiModel: Deterministic Compartmental Models with epiDCM'),  
  
  # Sidebar
  sidebarPanel(
  
    h4('Model Type'),
    selectInput(inputId='modtype', label='', 
                choices=c('SI', 'SIR', 'SIS')),
    br(), br(),
    
    h4('Initial Conditions'),
    numericInput(inputId='s.num', label='Number Susceptible', 
                 value=1000, min=0),
    numericInput(inputId='i.num', label='Number Infected', 
                 value=1, min=0),
    numericInput(inputId='r.num', label='Number Recovered', 
                 value=0, min=0),
    br(), br(),
    
    h4('Time'),
    numericInput(inputId='nsteps', label='Timesteps', 
                 value=500, min=0),
    numericInput(inputId='dt', label='dt (TS size)',
                 value=1, min=0),
    br(), br(),
    
    h4('Parameters'),
    numericInput(inputId='trans.rate', label='Transmission per Contact',
                min=0, max=1, value=1),
    br(),
    numericInput(inputId='act.rate', label='Acts per Timestep',
                 min=0, value=1),
    br(),
    numericInput(inputId='rec.rate', label='Recovery Rate (1/duration)',
                 min=0, value=0),
    br(),
    
    checkboxInput(inputId='vital', label=strong('Vital Dynamics'),
                  value=FALSE),
    conditionalPanel(
      condition='input.vital == true',
       numericInput(inputId='b.rate', label='Birth Rate',
                   min=0, value=0.0),
       br(),
       numericInput(inputId='ds.rate', label='Death Rate (Sus.)',
                   min=0, value=0.0),
       br(),
       numericInput(inputId='di.rate', label='Death Rate (Inf.)',
                   min=0, value=0.0),
       br(),
       numericInput(inputId='dr.rate', label='Death Rate (Rec.)',
                   min=0, value=0.0)
    )
  ),
  
  # Main panel
  mainPanel(
    tabsetPanel(
      tabPanel('Plot', 
               h4('Plot of Model Results'),
               plotOutput(outputId='MainPlot'),
               br(),
               selectInput(inputId='compsel', label='Plot Selection', 
                           choices=c('State Prevalence', 'State Size',
                                     'Incidence')),
               br(),
               sliderInput(inputId='alpha', label='Line Transparency',
                           min=0.1, max=1, value=0.8, step=0.1),
               br(),
               downloadButton('dlMainPlot', 'Download Plot')            
      ),
      tabPanel('Summary',
               h4('Summary from the Model'),
               numericInput(inputId='summTs', label='Timestep', value=1, min=1, max=500),
               numericInput('summDig', label='Significant Digits', value=3, min=0, max=8),
               verbatimTextOutput('outSummary'),
               br(),br(),
               
               h4('Compartment Plot'),
               plotOutput(outputId='CompPlot'),
               downloadButton('dlCompPlot', 'Download Plot'),
               br(),br()
      ),
      tabPanel('Data', 
               h4('Model Data'),
               downloadButton('dlData', 'Download Data'),
               br(),br(),
               tableOutput('outData')
      ),
      tabPanel('About', 
               p('This application solves and plots a deterministic, compartmental 
                  epidemic models. Models here are limited to one-group models, but 
                 two-group models are available directly through the epiDCM function. 
                 The underlying modeling software for this application is the', 
                 a('EpiModel', href='http://cran.r-project.org/web/packages/EpiModel/index.html'), 
                  'package in R. The web application is built with', 
                 a("Shiny.", href="http://www.rstudio.com/shiny/")),
               br(),
               strong('Author'), p('Samuel M. Jenness, Department of Epidemiology, 
                                    University of Washington')
      )
    )
  )
))

