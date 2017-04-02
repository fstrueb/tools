library(shiny)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "some cReative name"),
  dashboardSidebar(sidebarMenu(
    menuItem(
      text = 'Regulatory Element Scan',
      icon = icon('magic'),
      menuSubItem(text = 'Define custom range',
        tabName = 'defineRange'),
      menuSubItem(text = 'RE prediction from CAGE data',
        tabName = 'rePred'),
      menuSubItem(text = 'Scan Summary & Details',
        tabName = 'reSummary')
    ),
    menuItem(
      text = 'Regulatory Element Plot',
      tabName = 'rePlot',
      icon = icon('area-chart')
    ),
    menuItem(
      text = 'Regulatory Element Interactions',
      tabName = 'reIntx',
      icon = icon('exchange')
    )
  )),
  dashboardBody(tabItems(
    tabItem(
      ############# DEFINE RANGE ###########
      tabName = 'defineRange',
      fluidRow(
        box(
          title = strong('Define region of interest'),
          hr(),
          selectInput(
            inputId = 'chromSel',
            label = h5('Chromosome'),
            choices = list(
              'Chr1' = 'chr1',
              'Chr2' = 'chr2',
              'Chr3' = 'chr3',
              'Chr4' = 'chr4',
              'Chr12' = 'chr12'
            ),
            selected = 'chr2'
          ),
          ########## implement function that gives list according to chosen assembly/species
          numericInput(
            inputId = 'from',
            label = h5('From'),
            value = 109886546,
            min = 0,
            ######### define limits function for each chromosome later on
            max = 1e8
          ),
          numericInput(
            inputId = 'to',
            label = h5('To'),
            value = 109892833,
            min = 0,
            max = 1e8
          ),
          hr(),
          shinyjs::useShinyjs(),
          actionButton(inputId = 'acceptRange',
            label = 'Submit'),
          width = 3
        ),
        ####### Set the scanning parameters ########
        box(
          checkboxGroupInput(
            inputId = 'rangeSpecies',
            label = 'Select species',
            choices = c(
              'M. musculus' = '10090',
              'H. sapiens' = '9606',
              'A. thaliana' = '3702',
              'S. cerevisiae' = '4932',
              'D. melanogaster' = '7227',
              'C. elegans' = '6239',
              'R. norvegicus' = '10116'
            ),
            selected = 'M. musculus'
          ),
          checkboxGroupInput(
            inputId = 'rangeCollection',
            label = 'JASPAR Collection',
            choices = c(
              'CORE',
              'CNE',
              'PHYLOFACTS',
              'SPLICE',
              'POLII',
              'FAM',
              'PBM',
              'PBM_HOMEO',
              'PBM_HLH'
            ),
            selected = 'CORE'
          ),
          width = 3
        ),
        box(
          ##### dynamically update the motif choices #########
          htmlOutput('foundMotifsUI'),
          width = 4
        )
        ),
        fluidRow(box(
          dataTableOutput(outputId = 'rangeSummary'),
          width = 12
        ))),
        tabItem(tabName = 'rePred',
          fluidRow(box(title = 'predict'))),
        tabItem(tabName = 'reSummary',
          fluidRow(box(
            dataTableOutput(outputId = 'reSummaryTable'),
            width = 12
          ))),
        tabItem(tabName = 'rePlot',
          verticalLayout(
            box(
              title = strong('Select motifs to plot'),
              hr(),
              # selectInput(inputId = 'in1',
              #             label = 'choose motif',
              #             choices = list(),
              #             multiple = T,
              #             selectize = T),
              htmlOutput('selectUI'),
              actionButton(inputId = 'action1', label = 'Plot!'),
              plotOutput("results", height = 250),
              width = 12
            ),
            box(title = strong('TFBS results'),
              hr(),
              width = 10)
          ))
      ))
  )

    