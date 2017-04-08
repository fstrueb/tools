library(shiny)
library(dplyr)
library(rtracklayer)
library(DT)
library(shinyjs)

path = getwd()
source('../R_functions/plotPromoter.R')
source('../R_functions/makeRangeOfInterest.R')
source('../R_functions/scanRangeForTFBS.R')
source('../R_functions/siteSetToDataFrame.R')
source('../R_functions/unlistJASPAR.R')
source('../R_functions/standardizeSeqlevels.R')
source('../R_functions/siteSetListSummary.R')


shinyServer(function(input, output, session) {
  
  ####### DATA CHOICE #######
  reactive({
    if (input$rangeChoice == 'upload') {
      shinyjs::disable('chromSel')
      shinyjs::disable('from')
      shinyjs::disable('to')
    }
  })
  ############ FILE IMPORT #############
  rootVolumes = c(Home = normalizePath("~"), getVolumes()(), WD = '.')
  shinyFileChoose(input, 'file', roots = rootVolumes, session = session)
  loadCSVObj <- reactive({
    loadCSVFile(req(as.character(
      parseFilePaths(rootVolumes,
        input$file)$datapath
    )), input$csvsample1, input$csvsample2)
    
  })
  
  ############ REACTIVE VALUES ###########
  
  
  rangeObj = reactiveValues(
    # species = 'mouse' 
    # assembly = 'mm10'
    # standard: Lin7c-Ccdc34 locus, mm10
    assembly = 'mm10',
    chromosome = 'chr2',
    start = 109886546,
    end = 109892833)
  
  scanDetails = reactiveValues(
    species = 'mouse',
    collection = 'CORE')
  
  scanMotifs = reactiveValues(
    motif_ID = c('MA0001.1', 'MA0002.1'))
  
  scanResults = reactiveValues(
    df = data.frame(),
    summary = data.frame())
  
  
  ############ MESSAGE MENU ##############
  output$notificationMenu = renderMenu({
    notifications = list()
    if (is.null(scanResults$df)) {
      notifications[[1]] = notificationItem(text = 'Range not present, please provide one.', icon = icon('warning'))
    } else {
      notifications[[1]] = notificationItem(text = 'No notifications.', icon = icon('check'))
    }
    dropdownMenu(type = 'notification', .list = notifications)
  })
  
  ####### ISOLATE ###########
  # isolate checkbox button inputs 
  isolate({
    #rangeObj$assembly = input$annoChoice
    scanDetails$species = input$rangeSpecies
    scanDetails$collection = input$rangeCollection
    scanMotifs$motif_ID = input$resultMotifs
  })
  
  ########### SUBMIT AND SCAN ##########
  observeEvent(input$acceptRange, {
    # disable the button once the submit button is clicked
    shinyjs::disable('acceptRange')
    output$rangeSummary = renderDataTable(options = list(pageLength = 25, scrollX = T), {
      #input$acceptRange 
      # Create a Progress object
      progress <- shiny::Progress$new(style = 'notification')
      progress$set(message = "Scanning range...", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 4
        }
        progress$set(value = value, detail = detail)
      }
      isolate({
        rangeObj$assembly = input$annoChoice
        rangeObj$chromosome = input$chromSel
        rangeObj$start = input$from
        rangeObj$end = input$to
      })
      # scanMotifs$motif_ID = motif.list$motif_ID
      range = makeRange(chromosome = rangeObj$chromosome, from = rangeObj$start, to = rangeObj$end)
      cat('worked: ', as.character(GenomicRanges::seqnames(range)), '\n')
      ######### MAIN FUNCTIONS FOR PROGRESS BAR ##########
      withProgress(message = 'Scanning...', style = 'notification', value = 0.1, {
        ###### check if only specific JASPAR IDs are selected 
        isolate({
          if (is.null(input$resultMotifs)) {
            motif.list = unlistJASPAR(species = scanDetails$species, collection = scanDetails$collection)
            scanRes = scanRangeForTFBS(query = range, 
              motif.list = motif.list$motif_ID, 
              input.assembly = rangeObj$assembly, 
              updateProgress)
          } else {
            scanRes = scanRangeForTFBS(query = range, 
              motif.list = isolate(input$resultMotifs), 
              input.assembly = range$assembly, 
              updateProgress)
          }
        })
      })
      scanResSum = siteSetListSummary(query = range, siteSetList = scanRes)
      
      
      ### implement pvalue and motif return functionality here, update second progress bar
      
      # re-enable the submit button
      shinyjs::enable('acceptRange')
      # update the scanRes object
      scanResults$summary = scanResSum
      scanResSum
    }
      # enable scroll bar for data table output
    )
    # must update reSummary output
    output$reSummaryTable = renderDataTable(scanResults$df, options = list(scrollX = TRUE))
  })
  
  ############# dynamically change motif ID to select from for scan ##########
  foundMotifs = reactive({
    scanDetails$species = input$rangeSpecies
    scanDetails$collection = input$rangeCollection
    motif.list = unlistJASPAR(species = scanDetails$species, collection = scanDetails$collection)
    # if (is.null(motif.list)) {
    #   warning('motif.list is empty')
    # } else {
    #   motif.list = 'nothing defined'
    # }
    scanMotifs$motif_ID = motif.list
    motif.list$motif_ID
  })
  output$foundMotifsUI = renderUI({
    selectizeInput('resultMotifs', "Select JASPAR IDs", foundMotifs(), multiple = T)
  })
  
  ######## dynamically change motif ID input boxes in plot ##############
  searchResult = reactive({
    #scanResults$df %>% group_by(motif_ID) %>% tally() %>% arrange(-n) %>% mutate(n = as.character(n)) %>% tidyr::unite(col = result, everything(), sep = ', hits = ')
    unique(scanResults$df$motif_ID)
  })
  output$selectUI = renderUI({ 
    selectizeInput("plotMotifChoices", "Select your choice", searchResult())
  })
  
  ######### PLOT PART #############
  observeEvent(input$plotButton, {
    output$plotResults = renderPlot({
      if (is.null(scanResults$df)) {
        print('bla')
        plot(mtcars)
      } else {
        plotPromoter(
          range = scanResults$df,
          chr = 'chr12',
          motif.selection = isolate(input$plotMotifChoices)
        )
        plot
      }
    }) 
  })
  
  
})
