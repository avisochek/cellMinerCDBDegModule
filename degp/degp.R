source("degp/global.R", local = TRUE)
source("degp/ui.R", local = TRUE)
source("degp/analysis_outputs.R", local = TRUE)

### Server --------------------------------------------------------------------------------------------------
degpServer <- function(id, srcContentReactive) {
  moduleServer(id, function(input, output, session) {
  
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
    # Reactive to store limma fit
    degFit = reactiveVal(NULL),
    # Reactive to store gene ranks for fgsea
    fgseaStats = reactiveValues(geneRanking = NULL)
  )

  #Initialize data with Group column 
  observe({
    state$degFit(NULL)
    state$sampleData(initializeSampleData(srcContentReactive()[[input$dataSet]][["sampleData"]]))
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
  observeEvent(input$addGroup, {
    selectedRows <- input$dataSetTable_rows_selected
    groupName <- stringr::str_trim(input$groupName)
    
    if(length(selectedRows) > 2 && groupName != ""){
      
      #Update Group column 
      sampleTable <- state$sampleData()
      # Maximum of 3 allowed groups
      if(length(unique(na.omit(sampleTable$Group))) > 2) {
        shiny::showNotification("You may add a maximum of 3 groups. Please delete an existing group before creating another group", type = "error")
        return()
      }

      sampleTable$Group[selectedRows] <- groupName
      state$sampleData(sampleTable)
      state$degFit(NULL)
      
      groupCount <- sum(sampleTable$Group == groupName, na.rm = TRUE)
      
      #Reset selection in table
      DT::selectRows(dataSetTable_proxy, NULL)
      
      #Group creation notification
      shiny::showNotification(paste("Group", groupName, "now has", groupCount, "cell lines"))
      
      #Clear UI
      updateSelectizeInput(session, "selectIn1", selected = "")
      updateSelectizeInput(session, "selectIn2", selected = "")
      updateTextInput(session, "groupName", value = "")
      updateSelectizeInput(session, "tissueGroup", selected = character(0), choices = NULL)
    }
    else {
      shiny::showNotification("Please select at least 3 cell lines to add a group. If you used the tissue selector, you may need to add additional tissue groups or choose a tissue group that has more cell lines.", type = "error")
    }
    
  })
  
  observeEvent(input$deleteGroup, {
    groupName <- input$deleteGroup
    sampleTable <- state$sampleData()
    sampleTable$Group[sampleTable$Group == groupName] <- NA
    state$sampleData(sampleTable)
    state$degFit(NULL)
    shiny::showNotification(paste("Group", groupName, "deleted"))
  })
  
  output$groupInfoDisplay <- renderUI({
    groupCounts <- table(na.omit(state$sampleData()$Group))

    if (length(groupCounts) > 0) {
      do.call(tagList, lapply(names(groupCounts), function(gname) {
        fluidRow(
          column(5, strong(gname)),
          column(5, paste(groupCounts[[gname]], "cell lines")),
          column(2, tags$button(
            type = "button",
            class = "btn btn-link btn-xs",
            title = paste("Delete", gname),
            onclick = sprintf(
              "Shiny.setInputValue('%s', %s, {priority: 'event'})",
              session$ns("deleteGroup"),
              jsonlite::toJSON(gname, auto_unbox = TRUE)
            ),
            "x"
          ))
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
    selectizeInput(session$ns("tissueGroup"), "Select Tissue", choices = tissueSelectionOptions(), selected = character(0), multiple = TRUE)
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
  
  observe({
    groupChoices <- c("", getUniqueGroups())
    updateSelectizeInput(session, "selectIn1", choices = groupChoices)
    updateSelectizeInput(session, "selectIn2", choices = groupChoices)
  })
  
  
  # Run the analysis and return a datatable of results
  observeEvent(input$runAnalysis,{
    if (!is.null(input$selectIn1) && input$selectIn1 != "" && input$selectIn1 == input$selectIn2) {
      shiny::showNotification("Please select two different groups to run contrast", type = "error")
      return()
    }
    
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
        molPharmData <- srcContentReactive()[[input$dataSet]][["molPharmData"]]
        exprData <- if (!is.null(molPharmData$xsq)) molPharmData$xsq else molPharmData$exp
        controlCellLines <- sampleTable$Name[which(sampleTable$Group == input$selectIn1)]
        testCellLines <- sampleTable$Name[which(sampleTable$Group == input$selectIn2)]
        exprData1 <- exprData[, controlCellLines, drop = FALSE]
        exprData2 <- exprData[, testCellLines, drop = FALSE]
        
        #Increment progress 
        incProgress(0.1, detail = "Calculating differential expression...")
        
        #Calculate differential expression
        if (is.null(state$degFit())) {
          groupedRows <- !is.na(sampleTable$Group)
          state$degFit(calculateLimmaFit(
            exprData[, sampleTable$Name[groupedRows], drop = FALSE],
            sampleTable$Group[groupedRows]
          ))
        }
        degResults <- calculateLimmaContrast(state$degFit(), input$selectIn1, input$selectIn2)
        # Order by fold change, then p-value
        degResults <- degResults[with(degResults, order(-logFC, P.Value)), ]
        
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
        geneRanking <- degResults$logFC / ((degResults$P.Value) + 1)
        names(geneRanking) <- sub("^(xsq|exp)", "", rownames(degResults))
        
        state$fgseaStats$geneRanking <- sort(geneRanking, decreasing = TRUE)
        #Perform FGSEA analysis
        fgseaResults <- fgsea(pathways = combinedGeneSets,
                              stats = state$fgseaStats$geneRanking,
                              scoreType = "std")
        
        fgseaResults <- fgseaResults[with(fgseaResults, order(-NES)), ]
        fgseaResults <- fgseaResults[, c("pathway", "pval", "padj", "ES", "NES", "leadingEdge")] # add leading edge, add datatable scroll
        
        renderAnalysisOutputs(input, output, degResults, exprData1, exprData2, fgseaResults)
        
        incProgress(0.2, detail = "Finalizing analysis...")
      }
      else {
        shiny::showNotification("Please select two groups to run contrast", type = "error")
      }
      
    })
    
  })
  
  # New selection
  observeEvent(input$newSelection, {
    #Reset expression data
    state$degFit(NULL)
    state$fgseaStats$geneRanking <- NULL
    
    # #Reset sample data to its initial state without selections
    # state$sampleData(initializeSampleData(srcContentReactive()[[input$dataSet]][["sampleData"]]))

    #Reset all UI elements 
    # updateSelectizeInput(session, "selectIn1", selected = "")
    # updateSelectizeInput(session, "selectIn2", selected = "")
    # updateTextInput(session, "groupName", value = "")
    # updateSelectizeInput(session, "tissueGroup", selected = character(0), choices = NULL)

    #Return to first tab 
    showTab("mainTabset", "Input", select = TRUE)
    hideTab("mainTabset", "Results")
    hideTab("mainTabset", "Heatmap")
    hideTab("mainTabset", "Volcano Plot")
    hideTab("mainTabset", "Pathway Analysis")
    
    #Clear tables and plots 
    output$resultsTable <- DT::renderDT({datatable(data.frame())})  
    output$volcanoPlot <- renderPlot({NULL})  
    output$pathwayAnalysisResults <- DT::renderDT({datatable(data.frame())}) 
    output$pathwayAnalysisTopDotPlot <- renderPlot({NULL}) 
    
    #Clear data table 
    DT::selectRows(dataSetTable_proxy, NULL)
    
  })
  
  
  
  })
}
