source("global.R", local = TRUE)
source("ui.R", local = TRUE)
source("analysis_outputs.R", local = TRUE)

### Server --------------------------------------------------------------------------------------------------
server <- function(input, output,session) {
  
  hideTab("mainTabset", "Results")
  hideTab("mainTabset", "Heatmap")
  hideTab("mainTabset", "Volcano Plot")
  hideTab("mainTabset", "Pathway Analysis")
  
  initializeSampleData <- function(sampleTable) {
    sampleTable$Group <- NA
    sampleTable
  }

  state <- list(
    #Reactive to store sample data with group assignments
    sampleData = reactiveVal(),
    # Reactive to store gene ranks for fgsea
    fgseaStats = reactiveValues(geneRanking = NULL)
  )
  
  #Initialize data with Group column 
  observe({
    state$sampleData(initializeSampleData(srcContent[[input$dataSet]][["sampleData"]]))
  })
  
  #Note to self: Filter datasets (p2) 
  output$dataSetTable <- DT::renderDT({
    datatable(state$sampleData(),
              #Column selection below to restrict to cell line, tissue type information (source + OncoTree), and Group.
              options = list(
                columnDefs=list(
                  list(visible = FALSE, targets = c(0,5,6,7:(ncol(state$sampleData())-1)))
                )
              ),
              filter = 'top',
              # Enable row selection
              selection = 'multiple')
    
  })
  
  #Create datatable proxy to manipulate the table
  dataSetTable_proxy <- DT::dataTableProxy('dataSetTable')
  
  #Cell line selection and group creation 
  observeEvent(input$createGroup, {
    selectedRows <- input$dataSetTable_rows_selected
    groupName <- input$groupName
    
    if(length(selectedRows) > 2 && groupName != ""){
      
      #Update Group column 
      sampleTable <- state$sampleData()
      sampleTable$Group[selectedRows] <- groupName
      state$sampleData(sampleTable)
      
      groupCount <- sum(sampleTable$Group == groupName, na.rm = TRUE)
      
      #Reset selection in table
      DT::selectRows(dataSetTable_proxy, NULL)
      
      #Group creation notification
      shiny::showNotification(paste("Group", groupName, "now has", groupCount, "cell lines"))
      
      #Clear UI
      updateSelectInput(session, "selectIn1", selected = "")
      updateSelectInput(session, "selectIn2", selected = "")
      updateTextInput(session, "groupName", value = "")
      updateSelectizeInput(session, "tissueGroup", selected = character(0), choices = NULL)
    }
    else {
      shiny::showNotification("Please select at least 3 cell lines to create a group", type = "error")
    }
    
  })
  
  output$groupInfoDisplay <- renderUI({
    groupCounts <- table(na.omit(state$sampleData()$Group))

    if (length(groupCounts) > 0) {
      do.call(fluidRow, lapply(names(groupCounts), function(gname) {
        tagList(
          column(6, strong(gname)),
          column(6, paste(groupCounts[[gname]], "cell lines"))
        )
      }))
    } else {
      tags$p("No groups have been created yet.")
    }
  })
  
  # Search by tissue type
  tissueSelectionOptions <- reactive({
    data <- state$sampleData()
    #choices <- unique(data$OncoTree1)
    uniqueTissues <- unique(paste(data$OncoTree1, data$OncoTree2, sep = ": "))
    uniqueOncoTree1 <- unique(data$OncoTree1)
    
    combinedChoices <- unique(c(uniqueOncoTree1, uniqueTissues))
    combinedChoices
  })
  
  output$tissueSelector <- renderUI({
    selectizeInput("tissueGroup", "Select Tissue", choices = tissueSelectionOptions(), selected = character(0), multiple = TRUE)
  })
  
  observeEvent(input$tissueGroup, {
    selectedTissues <- input$tissueGroup
    data <- state$sampleData()
    selectedRows <- integer(0)
    for (selection in selectedTissues) {
      if (selection %in% data$OncoTree1) {
        # Select based on OncoTree1
        selectedRows <- c(selectedRows, which(data$OncoTree1 == selection))
      } else {
        # Select based on OncoTree2 (assumes selection is in "OncoTree1: OncoTree2" format)
        oncoTree1 <- strsplit(selection, ": ")[[1]][1]
        oncoTree2 <- strsplit(selection, ": ")[[1]][2]
        selectedRows <- c(selectedRows, which(data$OncoTree1 == oncoTree1 & data$OncoTree2 == oncoTree2))
      }
    }
    # Update row selection in the DataTable
    DT::selectRows(dataSetTable_proxy, unique(selectedRows))
  })
  
  #Reactive helper function for getting user-defined group names
  
  #runs analysis but no counts
  getUniqueGroups <- reactive({
    data <- state$sampleData()
    if("Group" %in% colnames(data)) {
      uniqueGroups <- unique(data$Group)
      return(na.omit(uniqueGroups))
    }
    else {
      return(character(0))
    }
  })
  
  # Allow user to select group 1 for contrast
  observe({
    # isolate({
    output$choice1 <-renderUI({
      selectizeInput("selectIn1", "Select Control Group:", choices = c(" ", getUniqueGroups()))
    })
    # })
  })
  # Allow user to select group 2 for contrast
  observe({
    # isolate({
    output$choice2 <- renderUI({
      selectizeInput("selectIn2", "Select Test Group:", choices = c(" ", getUniqueGroups()))
    })
    # })
  })
  
  
  # Run the analysis and return a datatable of results
  observeEvent(input$runAnalysis,{
    hideTab("mainTabset", "Input")
    showTab("mainTabset", "Results", select = TRUE)
    showTab("mainTabset", "Heatmap")
    showTab("mainTabset", "Volcano Plot")
    showTab("mainTabset", "Pathway Analysis")
    
    #Add progress bar 
    
    withProgress(message = 'Analysis in progress', {
      if (!is.null(input$selectIn1) && !is.null(input$selectIn2)) {
        # Retrieve RNA-Seq data for groups
        sampleTable <- state$sampleData()
        rnaSeqData <- srcContent[[input$dataSet]][["molPharmData"]][["xsq"]]
        controlCellLines <- sampleTable$Name[which(sampleTable$Group == input$selectIn1)]
        testCellLines <- sampleTable$Name[which(sampleTable$Group == input$selectIn2)]
        rnaSeqData1 <- rnaSeqData[, controlCellLines, drop = FALSE]
        rnaSeqData2 <- rnaSeqData[, testCellLines, drop = FALSE]
        
        #Increment progress 
        incProgress(0.1, detail = "Calculating differential expression...")
        
        #Calculate differential expression
        degResults <- calculateDEG(rnaSeqData1, rnaSeqData2)
        # Order by p_value, then fold change
        degResults <- degResults[with(degResults, order(P_Value, -Log2_Fold_Change)), ]
        
        # #Increment progress
        incProgress(0.2, detail = "Rendering results...")
        
        output$downloadResults <- downloadHandler(
          filename = function() {
            paste("DEG_results_", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            write.csv(degResults, file, row.names = FALSE)
          }
        )
        
        
        # FGSEA analysis
        incProgress(0.5, detail = "Performing pathway analysis...")
        
        #selectedGeneSet <- switch(input$geneSetChoice,
                                  #"hallmark" = hallmarkGeneSets,
                                  #"dtb" = dtbGeneSets)
        
        
        #Calculate FGSEA ranking metric as fold change over (p-value + 1)
        geneRanking <- degResults$Log2_Fold_Change / ((degResults$P_Value) + 1)
        names(geneRanking) <- degResults$Gene
        
        state$fgseaStats$geneRanking <- sort(geneRanking, decreasing = TRUE)
        #Perform FGSEA analysis
        fgseaResults <- fgsea(pathways = combinedGeneSets,
                              stats = state$fgseaStats$geneRanking,
                              scoreType = "pos")
        
        fgseaResults <- fgseaResults[with(fgseaResults, order(-NES)), ]
        fgseaResults <- fgseaResults[, c("pathway", "pval", "padj", "ES", "NES", "leadingEdge")] # add leading edge, add datatable scroll
        
        renderAnalysisOutputs(input, output, degResults, rnaSeqData1, rnaSeqData2, fgseaResults)
        
        incProgress(0.2, detail = "Finalizing analysis...")
      }
      else {
        shiny::showNotification("Please select two groups to run contrast", type = "error")
      }
      
    })
    
  })
  
  # Reset selections
  observeEvent(input$resetSelections, {
    #Reset expression data
    state$fgseaStats$geneRanking <- NULL
    
    #Reset sample data to its initial state without selections
    state$sampleData(initializeSampleData(srcContent[[input$dataSet]][["sampleData"]]))
    
    #Reset all UI elements 
    updateSelectInput(session, "selectIn1", selected = "")
    updateSelectInput(session, "selectIn2", selected = "")
    updateTextInput(session, "groupName", value = "")
    updateSelectizeInput(session, "tissueGroup", selected = character(0), choices = NULL)
    
    #Return to first tab 
    #updateTabsetPanel(session, "mainTabset", selected = "Input")
    showTab("mainTabset", "Input")
    hideTab("mainTabset", "Results")
    # hideTab("mainTabset", "Heatmap")
    hideTab("mainTabset", "Volcano Plot")
    hideTab("mainTabset", "Pathway Analysis")
    
    #Clear tables and plots 
    output$resultsTable <- DT::renderDT({datatable(data.frame())})  
    output$volcanoPlot <- renderPlot({NULL})  
    output$pathwayAnalysisResults <- DT::renderDT({datatable(data.frame())}) 
    
    #Clear data table 
    DT::selectRows(dataSetTable_proxy, NULL)
    
  })
  
  
  
}
shinyApp(ui, server)
