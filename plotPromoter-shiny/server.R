library(shiny)
library(dplyr)
library(rtracklayer)
library(DT)
library(shinyjs)
source('../plotPromoter.R')
source('../makeRangeOfInterest.R')
source('../scanRangeForTFBS.R')
source('../unlistJASPAR.R')

shinyServer(function(input, output, session) {
  
  ############ Define reactive values ###########
  rangeObj = reactiveValues(
    # species = 'mouse' 
    # assembly = 'mm10'
    # standard: Lin7c-Ccdc34 locus, mm10
    chromosome = 'chr2',
    start = 109886546,
    end = 109892833)
  
  scanDetails = reactiveValues(
    species = 'mouse',
    collection = 'CORE')
  
  scanMotifs = reactiveValues(
    motif_ID = c('MA0001.1', 'MA0002.1'))
  
  scanResults = reactiveValues(
    df = data.frame())
  
  # isolate checkbox button inputs 
  isolate({
    rangeObj$chromosome= input$chromSel
    rangeObj$start = input$from
    rangeObj$end = input$to
    scanDetails$species = input$rangeSpecies
    scanDetails$collection = input$rangeCollection
    scanMotifs$motif_ID = input$resultMotifs
  })
  
  ########### Submit and scan ##########
  observeEvent(input$acceptRange, {
    # disable the button once the submit button is clicked
    shinyjs::disable('acceptRange')
    output$rangeSummary = renderDataTable({
      # input$acceptRange 
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
      
      
      # scanMotifs$motif_ID = motif.list$motif_ID
      range = makeRangeOfInterest(chromosome = rangeObj$chromosome, from = rangeObj$start, to = rangeObj$end)
      
      ######### MAIN FUNCTIONS FOR PROGRESS BAR ##########
      withProgress(message = 'Scanning...', style = 'notification', value = 0.1, {
        ###### check if only specific JASPAR IDs are selected 
        isolate({
        if(is.null(input$resultMotifs)) {
          motif.list = unlistJASPAR(species = scanDetails$species, collection = scanDetails$collection)
          scanRes = scanRangeForTFBS(query = range, 
            motif.list = motif.list$motif_ID, 
            input.assembly = 'mm10', 
            return.p.val = F, 
            return.sequence = F, 
            updateProgress)
        } else {
          scanRes = scanRangeForTFBS(query = range, 
            motif.list = isolate(input$resultMotifs), 
            input.assembly = 'mm10', 
            return.p.val = F, 
            return.sequence = F, 
            updateProgress)
        }
      })
      })
      # re-enable the submit button
      shinyjs::enable('acceptRange')
      # update the scanRes object
      scanResults$df = scanRes
      scanRes
    },
      # enable scroll bar for data table output
      options = list(scrollX = TRUE))
    # must update reSummary output
    output$reSummaryTable = renderDataTable(scanResults$df, options = list(scrollX = TRUE))
  })
  
  ############# dynamically change motif ID to select from for scan ##########
  foundMotifs = reactive({
    scanDetails$species = input$rangeSpecies
    scanDetails$collection = input$rangeCollection
    motif.list = unlistJASPAR(species = scanDetails$species, collection = scanDetails$collection)
    scanMotifs$motif_ID = motif.list
    motif.list$motif_ID
  })
  output$foundMotifsUI = renderUI({
    selectizeInput('resultMotifs', "Select JASPAR IDs", foundMotifs(), multiple = T)
  })
  
  ############# dynamically change motif ID input boxes in plot ##############
  searchResult<- reactive({
    scanResults$df %>% group_by(motif_ID) %>% tally() %>% arrange(-n) %>% mutate(n = as.character(n)) %>% tidyr::unite(col = result, everything(), sep = ', hits = ')
  })
  output$selectUI <- renderUI({ 
    selectizeInput("result", "Select your choice", searchResult())
  })
  
  
  
  
})
