library(shiny)
library(shinyFiles)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "some cReative name", dropdownMenuOutput(outputId = 'notificationMenu')),
  dashboardSidebar(sidebarMenu(
    menuItem(
      text = 'Regulatory Element Scan',
      icon = icon('magic'),
      menuSubItem(text = 'Define custom range',
        tabName = 'defineRange'),
      menuSubItem(text  = 'Upload a set of ranges',
        tabName = 'uploadRange'),
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
        column(width = 6, 
          box(
            h5('which annotation?'),
            selectInput(
              inputId = 'annoChoice',
              label = 'Choose annotation:',
              choices = c(
                'mm9',
                'mm10',
                'hg37',
                'hg38'
              )
            ),
            hr(),
            h5('upload or define?'),
            radioButtons(
              inputId = 'rangeChoice',
              label = 'Choose method',
              choices = c(
                'Upload .csv file' = 'upload',
                'Define range manually' = 'manual'),
              selected = 'manual'),
            width = NULL
          ),
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
              selected = '10090'
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
            ), width = NULL
          )
        ),
        column(width = 6,
          conditionalPanel( 
            condition = "input.rangeChoice == 'manual'",
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
              width = NULL
            )),
          conditionalPanel(condition = "input.rangeChoice == 'upload'",
            box(
              shinyFilesButton('file', label = 'Select file', 'Provide a csv file containing ranges of interest', FALSE),
              width = NULL
            )),
          # box(
          #htmlOutput(outputId = 'rangeDefinitionUI'),
          # width = NULL,
          # title = strong('Upload ranges of interest'),
          # hr(),
          
          # ),
          ####### Set the scanning parameters ########
          # box(
          #   checkboxGroupInput(
          #     inputId = 'rangeSpecies',
          #     label = 'Select species',
          #     choices = c(
          #       'M. musculus' = '10090',
          #       'H. sapiens' = '9606',
          #       'A. thaliana' = '3702',
          #       'S. cerevisiae' = '4932',
          #       'D. melanogaster' = '7227',
          #       'C. elegans' = '6239',
          #       'R. norvegicus' = '10116'
          #     ),
          #     selected = 'M. musculus'
          #   ),
          #   checkboxGroupInput(
          #     inputId = 'rangeCollection',
          #     label = 'JASPAR Collection',
          #     choices = c(
          #       'CORE',
          #       'CNE',
          #       'PHYLOFACTS',
          #       'SPLICE',
          #       'POLII',
          #       'FAM',
          #       'PBM',
          #       'PBM_HOMEO',
          #       'PBM_HLH'
          #     ),
          #     selected = 'CORE'
          #   ),
          #   width = 3
          # ),
          box(
            ##### dynamically update the motif choices #########
            htmlOutput('foundMotifsUI'),
            width = NULL,
            hr(),
            shinyjs::useShinyjs(),
            actionButton(inputId = 'acceptRange',
              label = 'Submit'))
        ),
        fluidRow(box(
          dataTableOutput(outputId = 'rangeSummary'),
          width = 12
        )))),
    tabItem(tabName = 'rePred',
      fluidRow(box(title = 'predict'))),
    tabItem(tabName = 'reSummary',
      fluidRow(box(
        dataTableOutput(outputId = 'reSummaryTable'),
        width = 12
      ))),
    tabItem(tabName = 'rePlot',
      fluidRow(
        box(
          title = strong('Select motifs to plot'),
          hr(),
          # selectInput(inputId = 'in1',
          #             label = 'choose motif',
          #             choices = list(),
          #             multiple = T,
          #             selectize = T),
          htmlOutput('selectUI'),
          actionButton(inputId = 'plotButton', label = 'Plot!'),
          width = 3
        ),
        box(
          plotOutput(outputId = 'plotResults', height = 500),
          width = 8
        )
      )
    )
  ))
  
)

