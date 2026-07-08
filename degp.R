# Load packages ----------------------------------------------------------------
library(shiny)
library(rcellminer)
library(shinycssloaders)
library(plotly)
library(DT)
library(jsonlite)
library(stringr)
library(fgsea)
library(msigdbr)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(shinyalert)
library(shinyjs)
library(ComplexHeatmap)
library(heatmaply)
library(dplyr)

hasExpressionData <- function(dataSourceConfig) {
  packages <- dataSourceConfig[["packages"]]
  if (is.null(packages)) {
    return(FALSE)
  }

  any(vapply(packages, function(packageConfig) {
    molData <- packageConfig[["MolData"]]
    if (is.null(molData) || is.null(molData[["eSetListName"]])) {
      return(FALSE)
    }

    any(molData[["eSetListName"]] %in% c("xsq", "exp"), na.rm = TRUE)
  }, logical(1)))
}

filterAllNAExpressionColumns <- function(dataSourceContent) {
  for (assayName in c("xsq", "exp")) {
    expressionData <- dataSourceContent[["molPharmData"]][[assayName]]
    if (!is.null(expressionData)) {
      keepColumns <- colSums(is.na(expressionData)) != nrow(expressionData)
      dataSourceContent[["molPharmData"]][[assayName]] <- expressionData[, keepColumns, drop = FALSE]
    }
  }
  dataSourceContent
}

combineGeneSets <- function(set1, set2) {
  # Merge the two lists by combining gene sets with the same key and unioning genes in overlapping sets
  combinedSets <- lapply(union(names(set1), names(set2)), function(x) {
    unique(c(set1[[x]], set2[[x]]))
  })
  names(combinedSets) <- union(names(set1), names(set2))
  return(combinedSets)
}

# Format MSigDB gene sets for fgsea
formatMsigdbrGeneSets <- function(geneSetData) {
  geneSetData %>%
    distinct(gs_name, gene_symbol) %>%
    group_by(gs_name) %>%
    summarize(genes = list(gene_symbol), .groups = "drop") %>%
    tibble::deframe()
}

# Build Reactome pathway links from MSigDB metadata
buildReactomePathwayUrls <- function(reactomeData) {
  reactomeData %>%
    distinct(gs_name, gs_exact_source) %>%
    mutate(reactomeUrl = paste0("https://www.reactome.org/content/detail/", gs_exact_source)) %>%
    select(gs_name, reactomeUrl) %>%
    tibble::deframe()
}

# Define gene sets for pathway analysis 
hallmarkGeneSets <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
) %>%
  formatMsigdbrGeneSets()
dtbGeneSets <- gmtPathways("220411.Genelist.gmt")
reactomeGeneSetData <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
)
reactomeGeneSets <- formatMsigdbrGeneSets(reactomeGeneSetData)
reactomePathwayUrls <- buildReactomePathwayUrls(reactomeGeneSetData)

# Combine the gene sets
combinedGeneSets <- combineGeneSets(hallmarkGeneSets, dtbGeneSets)
combinedGeneSets <- combineGeneSets(combinedGeneSets, reactomeGeneSets)

### User Interface --------------------------------------------------------------------------------------------------
degpInput <- function(id) {
  ns <- NS(id)
  tabPanel("Differential Expression Analysis",
  fluidPage(
  useShinyjs(),
  useShinyalert(),
    tags$head(
      tags$style(HTML("
        :root {
          --app-message-font-size: 18px;
        }
        .shiny-notification {
          position: fixed !important;
          top: 40% !important;
          left: 50% !important;
          width: auto !important;
          padding: 10px !important;
        }
        .shiny-notification,
        .shiny-progress-notification,
        .shiny-progress-notification .progress-message,
        .shiny-progress-notification .progress-detail,
        .sweet-alert,
        .sweet-alert h2,
        .sweet-alert p,
        .swal-title,
        .swal-text {
          font-size: var(--app-message-font-size) !important;
        }
        .sidebar-step-heading {
          display: block;
          font-size: 22px;
          font-weight: 600;
          line-height: 1.25;
          margin: 0 0 10px;
        }
        .sidebar-helper-text {
          font-size: 14px;
          line-height: 1.35;
          margin: 0 0 6px;
        }
        .heatmap-output-container .shiny-output-error-validation {
          color: #A94442;
          font-size: 20px;
          font-weight: 500;
          line-height: 1.4;
        }
      "))
    ),
             sidebarLayout(
               #tags$head(tags$style(".shiny-notification {position: fixed; top: 60% ;left: 50%")),
               sidebarPanel(
                 selectInput(
                   ns("dataSet"),
                   label = tags$span(class = "sidebar-step-heading", "1. Select Dataset"),
                   choices = character(0)
                 ),
                 br(),
                 br(),
                 tags$div(class = "sidebar-step-heading", "2. Select Groups"),
                 tags$p(class = "sidebar-helper-text", "Type a name for the group."),
                 tags$p(class = "sidebar-helper-text", "Either select cell lines from the table, or choose 1 or more tissue types from the tissue selector."),
                 tags$p(class = "sidebar-helper-text", "When you are done, click 'Add Group'."),
                 br(),
                 textInput(ns("groupName"), "Group Name", value = ""),
                 fluidRow(
                   column(8, uiOutput(ns("tissueSelector"))),
                   column(
                     4,
                     tags$div(
                       style = "margin-top: 25px;",
                       actionButton(inputId = ns("addGroup"), label = "Add Group")
                     )
                   )
                 ),
                 br(),
                 br(),
                 uiOutput(ns("groupInfoDisplay")),
                 br(),
                 tags$div(
                   class = "sidebar-step-heading",
                   "3. Select Contrast",
                   actionLink(
                     ns("contrastHelp"),
                     label = NULL,
                     icon = icon("question-circle")
                   )
                 ),
                 selectizeInput(ns("selectIn2"), "Select Test Group:", choices = c("")),
                 selectizeInput(ns("selectIn1"), "Select Control Group:", choices = c("")),
                 br(),
                 tags$div(class = "sidebar-step-heading", "4. Run Analysis"),
                 actionButton(ns("runAnalysis"), "Run"),
                 actionButton(ns("newSelection"), "Change Selection")
               ),
               mainPanel(
                 tabsetPanel(id = ns("mainTabset"),
                             tabPanel("Input",
                                      fluidRow(
                                        DT::DTOutput(ns("dataSetTable"))
                                      )
                             ),
                             tabPanel("Results",
                                      uiOutput(ns("resultsHeading")),
                                      fluidRow(
                                        withSpinner(DT::DTOutput(ns("resultsTable"))),
                                        downloadButton(ns("downloadResults"), "Download Results")
                                      )),
                             tabPanel("Volcano Plot",
                                      uiOutput(ns("volcanoHeading")),
                                      plotOutput(ns("volcanoPlot")),
                                      fluidRow(
                                        column(
                                          4,
                                          sliderInput(
                                            ns("volcanoFdrThreshold"),
                                            "-Log10 FDR significance threshold",
                                            min = 1,
                                            max = 5,
                                            value = -log10(0.05),
                                            step = 0.1
                                          ),
                                          sliderInput(
                                            ns("volcanoLogFcThreshold"),
                                            "Log2 fold-change significance threshold",
                                            min = 0.1,
                                            max = 5,
                                            value = 1,
                                            step = 0.1
                                          ),
                                          actionButton(ns("updateVolcanoPlot"), "Update Plot")
                                        ),
                                        column(
                                          3,
                                          style = "margin-top: 25px;",
                                          uiOutput(ns("volcanoSignificanceCounts"))
                                        )
                                      )),
                             tabPanel("Heatmap",
                                      uiOutput(ns("heatmapHeading")),
                                      tags$div(
                                        class = "heatmap-output-container",
                                        plotlyOutput(
                                          ns("heatmapPlot"),
                                          width = "100%",
                                          height = "700px"
                                        )
                                      )),
                             tabPanel("Pathway Analysis",
                                      # fluidRow(
                                      #   column(12,
                                      #          div(style = "float: right;",
                                      #              selectInput("geneSetChoice", "Choose Gene Set for Pathway Analysis:",
                                      #                          choices = c("Hallmark Pathways Gene Set" = "hallmark",
                                      #                                      "DTB Gene Set" = "dtb"))
                                      #              )
                                      #          )
                                      # ),
                                      #actionButton("updatePathway", "Update Pathway Analysis"),
                                      tabsetPanel(
                                        tabPanel("FGSEA Table",
                                                 uiOutput(ns("pathwayAnalysisHeading")),
                                                 fluidRow(
                                                   DT::DTOutput(ns("pathwayAnalysisResults")),
                                                   downloadButton(ns("downloadFgseaResults"), "Download Results")
                                                 )),
                                        tabPanel("Top FGSEA Plot",
                                                 uiOutput(ns("pathwayAnalysisTopDotPlotHeading")),
                                                 plotOutput(ns("pathwayAnalysisTopDotPlot"), height = "700px"))
                                      ))
                 )
               )
             )
    )
  )
}

### Server --------------------------------------------------------------------------------------------------
degpServer <- function(input, output, session, srcContentReactive, config){
  
  hideTab("mainTabset", "Results")
  hideTab("mainTabset", "Heatmap")
  hideTab("mainTabset", "Volcano Plot")
  hideTab("mainTabset", "Pathway Analysis")
  
  initializeSampleData <- function(sampleTable) {
    sampleTable$Group <- NA
    sampleTable
  }

  filteredSrcContent <- reactiveVal(NULL)

  observeEvent(srcContentReactive(), {
    expressionConfig <- Filter(hasExpressionData, config)
    srcContent <- srcContentReactive()[names(expressionConfig)]

    srcContent <- lapply(srcContent, function(dataSource) {
      dataSource <- filterAllNAExpressionColumns(dataSource)
      molPharmData <- dataSource[["molPharmData"]]
      exprData <- if (!is.null(molPharmData$xsq)) molPharmData$xsq else molPharmData$exp
      dataSource[["sampleData"]] <- dataSource[["sampleData"]][dataSource[["sampleData"]]$Name %in% colnames(exprData), , drop = FALSE]
      dataSource
    })

    filteredSrcContent(srcContent)
    # Indicate in the dataset selection whether we are using rna-seq or microarray data for each dataset
    dataSourceChoices <- setNames(
      names(srcContent),
      vapply(names(srcContent), function(x) {
        mol_data <- srcContent[[x]][["molPharmData"]]
        assay_label <- if (!is.null(mol_data$xsq)) "RNA-seq" else "microarray"
        paste0(expressionConfig[[x]][["displayName"]], " (", assay_label, ")")
      }, character(1))
    )
    selectedDataSourceName <- if ("nci60" %in% dataSourceChoices) "nci60" else dataSourceChoices[[1]]
    updateSelectInput(session, "dataSet", choices = dataSourceChoices, selected = selectedDataSourceName)
  }, once = TRUE)

  selectedDataSource <- reactive({
    srcContent <- req(filteredSrcContent())
    req(input$dataSet)
    srcContent[[input$dataSet]]
  })

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
    state$sampleData(initializeSampleData(selectedDataSource()[["sampleData"]]))
  })
  
  output$dataSetTable <- DT::renderDT({
    sampleData <- req(state$sampleData())
    datatable(sampleData,
              #Column selection below to restrict to cell line, tissue type information (source + OncoTree), and Group.
              options = list(
                columnDefs=list(
                  list(visible = FALSE, targets = c(0,5,6,7:(ncol(sampleData)-1)))
                )
              ),
              filter = 'top',
              # Enable row selection
              selection = 'multiple')
    
  })

  observeEvent(input$contrastHelp, {
    shinyalert::shinyalert(
      title = "Contrast",
      text = paste(
        "Contrast defines the comparison used to calculate differential",
        "expression. It specifies which group is compared against another",
        "(e.g., test vs. control), determining the direction and",
        "interpretation of the reported statistics."
      ),
      type = "info"
    )
  })
  
  #Create datatable proxy to manipulate the table
  dataSetTable_proxy <- DT::dataTableProxy('dataSetTable')
  
  #Cell line selection and group creation 
  observeEvent(input$addGroup, {
    selectedRows <- input$dataSetTable_rows_selected
    groupName <- stringr::str_trim(input$groupName)
    
    if (groupName ==""){
      shinyalert::shinyalert(
        title = "Error",
        text = "Please type a name for the group.",
        type = "error"
      )
      return()
    }

    if (length(selectedRows) < 3) {
      shinyalert::shinyalert(
        title = "Error",
        text = paste(
          "Please select at least 3 cell lines to add a group.",
          "If you used the tissue selector, you may need to add additional",
          "tissue groups or choose a tissue group that has more cell lines."
        ),
        type = "error"
      )
      return()
    }
       
    sampleTable <- state$sampleData()
    if (groupName %in% na.omit(sampleTable$Group)) {
      shinyalert::shinyalert(
        title = "Error",
        text = paste(
          "A group named", groupName, "already exists.",
          "Please choose a different group name.",
          "If you would like to replace the existing group, please delete it first."
        ),
        type = "error"
      )
      return()
    }

    # Maximum of 3 allowed groups
    if(length(unique(na.omit(sampleTable$Group))) > 2) {
      shinyalert::shinyalert(
        title = "Error",
        text = paste(
          "You may add a maximum of 3 groups.",
          "Please delete an existing group before creating another group."
        ),
        type = "error"
      )
      return()
    }

    #Update Group column
    sampleTable$Group[selectedRows] <- groupName
    groupCounts <- table(na.omit(sampleTable$Group))
    undersizedGroups <- names(groupCounts[groupCounts < 3])
    if (length(undersizedGroups) > 0) {
      sampleTable$Group[sampleTable$Group %in% undersizedGroups] <- NA
    }
    state$sampleData(sampleTable)
    state$degFit(NULL)
    
    groupCount <- sum(sampleTable$Group == groupName, na.rm = TRUE)
    
    #Reset selection in table
    DT::selectRows(dataSetTable_proxy, NULL)
    
    #Group creation notification
    shiny::showNotification(paste("Group", groupName, "now has", groupCount, "cell lines"))
    if (length(undersizedGroups) > 0) {
      shiny::showNotification(
        paste(
          "Removed",
          paste(undersizedGroups, collapse = ", "),
          "because they had fewer than 3 cell lines."
        ),
        type = "warning"
      )
    }
    
    #Clear UI
    updateSelectizeInput(session, "selectIn1", selected = "")
    updateSelectizeInput(session, "selectIn2", selected = "")
    updateTextInput(session, "groupName", value = "")
    updateSelectizeInput(session, "tissueGroup", selected = character(0), choices = NULL)
    
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
    data <- req(state$sampleData())
    #choices <- unique(data$OncoTree1)
    uniqueTissues <- unique(paste(data$OncoTree1, data$OncoTree2, sep = ": "))
    uniqueOncoTree1 <- unique(data$OncoTree1)
    
    combinedChoices <- unique(c(uniqueOncoTree1, uniqueTissues))
    combinedChoices
  })
  
  output$tissueSelector <- renderUI({
    selectizeInput(session$ns("tissueGroup"), "Select Tissue(s)", choices = tissueSelectionOptions(), selected = character(0), multiple = TRUE)
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
      shinyalert::shinyalert(
        title = "Error",
        text = paste(
          "The test and control group can not be the same.",
          "Please select two different groups."
        ),
        type = "error"
      )
      return()
    }

    if ((is.null(input$selectIn1) | input$selectIn1 == "") |
        (is.null(input$selectIn2) | input$selectIn2 == "")) {
      shinyalert::shinyalert(
        title = "Error",
        text = "Please select two groups to run contrast.",
        type = "error"
      )
      return()
    }
    
    hideTab("mainTabset", "Input")
    showTab("mainTabset", "Results", select = TRUE)
    showTab("mainTabset", "Heatmap")
    showTab("mainTabset", "Volcano Plot")
    showTab("mainTabset", "Pathway Analysis")
    
    fgseaSeed <- 123
    
    #Add progress bar 
    
    withProgress(message = 'Analysis in progress', {
        # Retrieve RNA-Seq data for groups
        sampleTable <- state$sampleData()
        molPharmData <- selectedDataSource()[["molPharmData"]]
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
        degResults$`Mean expression test` <- rowMeans(
          exprData2[rownames(degResults), , drop = FALSE],
          na.rm = TRUE
        )
        degResults$`Mean expression ctrl` <- rowMeans(
          exprData1[rownames(degResults), , drop = FALSE],
          na.rm = TRUE
        )
        
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
        set.seed(fgseaSeed)
        fgseaResults <- fgsea(pathways = combinedGeneSets,
                              stats = state$fgseaStats$geneRanking,
                              nperm = 10000,
                              minSize = 5,
                              maxSize = 500,
                              scoreType = "std")
        
        fgseaResults <- fgseaResults[with(fgseaResults, order(-NES)), ]
        fgseaResults <- fgseaResults[, c("pathway", "pval", "padj", "ES", "NES", "leadingEdge")] # add leading edge, add datatable scroll

        output$downloadFgseaResults <- downloadHandler(
          filename = function() {
            paste("FGSEA_results_", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            fgseaDownloadResults <- as.data.frame(fgseaResults)
            fgseaDownloadResults$leadingEdge <- vapply(fgseaDownloadResults$leadingEdge, paste, character(1), collapse = ", ")
            write.csv(fgseaDownloadResults, file, row.names = FALSE)
          }
        )
        
        renderAnalysisOutputs(input, output, session, degResults, exprData1, exprData2, fgseaResults)
        
        incProgress(0.2, detail = "Finalizing analysis...")
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
    output$resultsHeading <- renderUI(NULL)
    output$volcanoPlot <- renderPlot({NULL})  
    output$volcanoHeading <- renderUI(NULL)
    output$volcanoSignificanceCounts <- renderUI(NULL)
    output$heatmapHeading <- renderUI(NULL)
    output$pathwayAnalysisResults <- DT::renderDT({datatable(data.frame())}) 
    output$pathwayAnalysisHeading <- renderUI(NULL)
    output$pathwayAnalysisTopDotPlot <- renderPlot({NULL}) 
    output$pathwayAnalysisTopDotPlotHeading <- renderUI(NULL)
    
    #Clear data table 
    DT::selectRows(dataSetTable_proxy, NULL)
    
  })
}  
  
  
### Analysis Outputs --------------------------------------------------------------------------------------------------
renderAnalysisOutputs <- function(input, output, session, degResults, exprData1, exprData2, fgseaResults) {
  output$resultsHeading <- renderUI({
    tagList(
      tags$h3(paste0("Differential Expression Analysis: ", input$selectIn2, " vs ", input$selectIn1)),
      helpText("Search table by regular expressions")
    )
  })

  # Render in results table
  output$resultsTable <- DT::renderDT({
    resultsColumnNames <- c(
      Gene = "Gene",
      logFC = "Log Fold Change",
      AveExpr = "Mean expression",
      t = "T Value",
      P.Value = "P Value",
      adj.P.Val = "Adjusted P Value",
      B = "B Statistic",
      Annotation = "Annotations",
      `Mean expression test` = "Mean expression test",
      `Mean expression ctrl` = "Mean expression ctrl"
    )
    hiddenColumns <- match(
      c("AveExpr", "t", "P.Value", "B"),
      colnames(degResults)
    ) - 1

    resultsTable <- datatable(degResults,
              rownames = FALSE,
              extensions = "Buttons",
              colnames = unname(resultsColumnNames[colnames(degResults)]),
              options = list(
                dom = "Bfrtip",
                buttons = list(
                  list(extend = "colvis", text = "Column visibility")
                ),
                search = list(regex = TRUE),
                columnDefs = list(
                  list(visible = FALSE, targets = hiddenColumns)
                )
              ),
              filter = 'top'
              
    )
    
    resultsTable <- DT::formatRound(
      resultsTable,
      columns = c(
        "logFC",
        "AveExpr",
        "t",
        "B",
        "Mean expression test",
        "Mean expression ctrl"
      ),
      digits = 2
    )
    resultsTable <- DT::formatSignif(resultsTable, columns = c("P.Value", "adj.P.Val"), digits = 3)
    resultsTable
  })
  
  
  output$volcanoHeading <- renderUI({
    tags$h3(paste0("Volcano Plot: ", input$selectIn2, " vs ", input$selectIn1))
  })

  volcanoResults <- degResults[abs(degResults$logFC) >= 0.10, ]
  volcanoLogFdrRange <- range(-log10(volcanoResults$adj.P.Val), na.rm = TRUE)
  volcanoLogFdrMax <- max(1, ceiling(volcanoLogFdrRange[2]))
  volcanoLogFcMax <- ceiling(max(abs(volcanoResults$logFC), na.rm = TRUE) * 10) / 10
  updateSliderInput(
    session,
    "volcanoFdrThreshold",
    min = 1,
    max = volcanoLogFdrMax,
    value = max(1, min(-log10(0.05), volcanoLogFdrMax))
  )
  updateSliderInput(
    session,
    "volcanoLogFcThreshold",
    min = 0.1,
    max = volcanoLogFcMax,
    value = min(1, volcanoLogFcMax)
  )

  volcanoThresholds <- eventReactive(input$updateVolcanoPlot, {
    fdrThresholdLog10 <- req(input$volcanoFdrThreshold)

    list(
      fdrThresholdLog10 = fdrThresholdLog10,
      fdrThreshold = 10 ^ -fdrThresholdLog10,
      logFcThreshold = req(input$volcanoLogFcThreshold)
    )
  }, ignoreNULL = FALSE)

  volcanoPlotData <- reactive({
    thresholds <- volcanoThresholds()
    fdrThreshold <- thresholds$fdrThreshold
    fdrThresholdLog10 <- thresholds$fdrThresholdLog10
    logFcThreshold <- thresholds$logFcThreshold
    fdrThresholdLog10Label <- format(signif(fdrThresholdLog10, 3), scientific = FALSE, trim = TRUE)
    logFcThresholdLabel <- format(signif(logFcThreshold, 3), scientific = FALSE, trim = TRUE)
    volcanoResults <- degResults[abs(degResults$logFC) >= 0.10, ]
    volcanoResults$P.Value[volcanoResults$adj.P.Val == 0] <- 10e-10

    log_FC = volcanoResults$logFC
    log_pval = -log10(volcanoResults$adj.P.Val)
    upregulatedLabel <- paste0(
      "Upregulated: Log2 FoldChange >= ",
      logFcThresholdLabel,
      ", -Log10 FDR >= ",
      fdrThresholdLog10Label
    )
    downregulatedLabel <- paste0(
      "Downregulated: Log2 FoldChange <= -",
      logFcThresholdLabel,
      ", -Log10 FDR >= ",
      fdrThresholdLog10Label
    )
    topGeneLabel <- "Top 10 upregulated / downregulated"

    isSignificant <- volcanoResults$adj.P.Val <= fdrThreshold
    Significant=rep("Not Significant",length(log_FC))
    Significant[which(isSignificant & volcanoResults$logFC>=logFcThreshold)]=upregulatedLabel
    Significant[which(isSignificant & volcanoResults$logFC<=-logFcThreshold)]=downregulatedLabel

    gene = sub("^xsq", "", rownames(volcanoResults))
    volcano_data=data.frame(gene,log_FC,log_pval,Significant)
    significantUpregulated <- volcano_data[volcano_data$Significant == upregulatedLabel, ]
    significantDownregulated <- volcano_data[volcano_data$Significant == downregulatedLabel, ]
    topUpregulatedGenes <- head(
      significantUpregulated[order(-significantUpregulated$log_FC, -significantUpregulated$log_pval), "gene"],
      10
    )
    topDownregulatedGenes <- head(
      significantDownregulated[order(significantDownregulated$log_FC, -significantDownregulated$log_pval), "gene"],
      10
    )
    volcano_data$TopGene <- ifelse(volcano_data$gene %in% c(topUpregulatedGenes, topDownregulatedGenes),
                                   topGeneLabel,
                                   NA)

    list(
      data = volcano_data,
      significantUpregulated = significantUpregulated,
      significantDownregulated = significantDownregulated,
      fdrThreshold = fdrThreshold,
      fdrThresholdLog10 = thresholds$fdrThresholdLog10,
      logFcThreshold = logFcThreshold,
      upregulatedLabel = upregulatedLabel,
      downregulatedLabel = downregulatedLabel,
      topGeneLabel = topGeneLabel
    )
  })

  output$volcanoPlot <- renderPlot({
    plotData <- volcanoPlotData()
    volcano_data <- plotData$data
    topGeneData <- volcano_data[!is.na(volcano_data$TopGene), ]
    pointSize <- 1.5
    xAxisMax <- max(2, abs(volcano_data$log_FC), na.rm = TRUE)
    yAxisMax <- max(2, volcano_data$log_pval, na.rm = TRUE) * 1.1
    volcanoColors <- c("deepskyblue", "red", "black", "grey")
    names(volcanoColors) <- c(
      plotData$downregulatedLabel,
      plotData$upregulatedLabel,
      plotData$topGeneLabel,
      "Not Significant"
    )
    volcano_plot <- ggplot(volcano_data, aes(x = log_FC, y = log_pval, color = Significant)) +
      geom_point(size = pointSize) +
      geom_point(
        data = topGeneData,
        aes(color = TopGene),
        shape = 1,
        fill = NA,
        size = pointSize,
        stroke = 0.8
      ) +
      ggrepel::geom_label_repel(
        data = topGeneData,
        aes(label = gene),
        color = "black",
        fill = "white",
        size = 3,
        label.size = 0.15,
        label.padding = grid::unit(0.1, "lines"),
        max.overlaps = Inf,
        seed = 123,
        show.legend = FALSE
      ) +
      geom_hline(yintercept = plotData$fdrThresholdLog10, linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-plotData$logFcThreshold, plotData$logFcThreshold), linetype = "dashed", color = "black") +
      scale_color_manual(
        name = "Significance",
        values = volcanoColors,
        breaks = names(volcanoColors)
      ) +
      labs(
        x = "Log2 Fold Change",
        y = "-Log10 FDR"
      ) +
      scale_x_continuous(
        limits = c(-xAxisMax, xAxisMax)
      ) +
      scale_y_continuous(
        limits = c(0, yAxisMax)
      ) +
      theme_classic() +
      theme(
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.key.height = grid::unit(1.2, "lines"),
        legend.key.width = grid::unit(0.8, "lines")
      )
    
    volcano_plot
    
  })

  output$volcanoSignificanceCounts <- renderUI({
    plotData <- volcanoPlotData()

    tags$p(
      tags$strong("Significant upregulated genes: "),
      nrow(plotData$significantUpregulated),
      tags$br(),
      tags$strong("Significant downregulated genes: "),
      nrow(plotData$significantDownregulated)
    )
  })
  
  
  output$heatmapHeading <- renderUI({
    tagList(
      tags$h3(paste0("Heatmap: ", input$selectIn2, " vs ", input$selectIn1)),
      helpText(
        "Top 10 significant upregulated and top 10 significant downregulated genes",
        tags$br(),
        "Distance metric: Pearson correlation",
        tags$br(),
        "Clustering method: complete linkage"
      )
    )
  })

  ## Heatmap
  heatmapPlotCache <- NULL

  output$heatmapPlot <- renderPlotly({
    if (!is.null(heatmapPlotCache)) {
      return(heatmapPlotCache)
    }
    
    significantUpregulated <- degResults[degResults$adj.P.Val < 0.05 & degResults$logFC > 1, ]
    significantDownregulated <- degResults[degResults$adj.P.Val < 0.05 & degResults$logFC < -1, ]

    shiny::validate(
      shiny::need(
        nrow(significantUpregulated) >= 3 && nrow(significantDownregulated) >= 3,
        "Not enough significant genes to display a heatmap. The heatmap requires at least three significant upregulated genes and three significant downregulated genes."
      )
    )

    withProgress(message = "Loading heatmap", {
      incProgress(0.1, detail = "Selecting significant genes...")

    topUpregulated <- head(significantUpregulated[order(-significantUpregulated$logFC, significantUpregulated$adj.P.Val), ], 10)
    topDownregulated <- head(significantDownregulated[order(significantDownregulated$logFC, significantDownregulated$adj.P.Val), ], 10)
    
    geneNames <- c(rownames(topDownregulated), rownames(topUpregulated))

      incProgress(0.2, detail = "Preparing expression data...")

    exprData1 <- data.frame(exprData1, check.names = FALSE) %>% select_if(~ !all(is.na(.)))
    exprData2 <- data.frame(exprData2, check.names = FALSE) %>% select_if(~ !all(is.na(.)))

    expressionData <- cbind(exprData1, exprData2)
    groupColors <- setNames(
      c("#0072B2", "#D55E00"),
      c(input$selectIn1, input$selectIn2)
    )
    columnGroups <- data.frame(
      Group = factor(
        c(
          rep(input$selectIn1, ncol(exprData1)),
          rep(input$selectIn2, ncol(exprData2))
        ),
        levels = names(groupColors)
      )
    )
    rownames(columnGroups) <- colnames(expressionData)
    
    expressionData <- expressionData[geneNames, , drop = FALSE]
    #Remove 'xsq' prefix from gene names
    rownames(expressionData) <- gsub("^xsq", "", rownames(expressionData))
    
    columnLabelSize <- round(min(16, max(8, 260 / ncol(expressionData))))
    # Height proportions for each component of the heatmap plot
    heatmap_body_height <- 0.8
    column_group_height <- 0.05
    dendrogram_height <- 0.15

      incProgress(0.5, detail = "Clustering heatmap...")

    heatmapPlot <- heatmaply(
      expressionData,
      xlab = "Cell Lines",
      ylab = "Genes",
      label_names = c("Gene", "Cell line", "Z score"),
      key.title = "Z score",
      colors = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(255),
      col_side_colors = columnGroups,
      col_side_palette = groupColors,
      dendrogram = "column",
      show_dendrogram = c(FALSE, TRUE),
      dend_hoverinfo = FALSE,
      distfun_col = "pearson",
      hclust_method = "complete",
      # Needed to make room for the legend
      margins = c(NA, NA, NA, 180),
      subplot_heights = c(dendrogram_height,column_group_height,heatmap_body_height),
      scale = "row",
      fontsize_col = columnLabelSize,
      fontsize_row = 10,
      custom_hovertext = matrix(
        paste0(
          "Original log2(FPKM+1): ",
          trimws(format(as.matrix(expressionData), digits = 4))
        ),
        nrow = nrow(expressionData),
        dimnames = dimnames(expressionData)
      )
    )

      incProgress(0.8, detail = "Finalizing heatmap...")

    ## Hide tooltip for group info display
    heatmapPlot$x$data <- lapply(heatmapPlot$x$data, function(trace) {
      if (identical(trace$yaxis, "y2")) {
        trace$hoverinfo <- "skip"
      }

      trace
    })
    heatmapPlot$x$layout$showlegend <- TRUE
    heatmapPlot$x$layout$legend$font <- list(size = 16)
    # Leave column group title blank
    heatmapPlot$x$layout$legend$title$text <- ""
    heatmapPlotCache <<- layout(
      heatmapPlot,
      autosize = TRUE
    )
    heatmapPlotCache
    })
    
  })
  
  
  # Display FGSEA results
  output$pathwayAnalysisHeading <- renderUI({
    tags$h3(paste0("Pathway Analysis Results: ", input$selectIn2, " vs ", input$selectIn1))
  })

  output$pathwayAnalysisResults <- renderDT({
   truncateFgseaCell <- function(value, link = NA_character_) {
     value <- paste(value, collapse = ", ")
     label <- htmltools::htmlEscape(value)
     if (!is.na(link) && nzchar(link)) {
       label <- sprintf("<a href='%s' target='_blank'>%s</a>", htmltools::htmlEscape(link, attribute = TRUE), label)
     }
     tooltip <- htmltools::htmlEscape(gsub("_", "_ ", value), attribute = TRUE)
     sprintf("<div style='max-width: 260px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;' data-toggle='tooltip' data-title='%s'>%s</div>", tooltip, label)
   }
   
   fgseaResults$pathway <- vapply(fgseaResults$pathway, function(pathway) {
     truncateFgseaCell(pathway, unname(reactomePathwayUrls[pathway]))
   }, character(1))
   fgseaResults$leadingEdge <- vapply(fgseaResults$leadingEdge, truncateFgseaCell, character(1))
   
   pathwayTable <- datatable(fgseaResults,
              escape = FALSE,
              callback = DT::JS("
                var activateFgseaTooltips = function() {
                  $(table.table().container()).find('[data-toggle=\"tooltip\"]').tooltip({
                    title: function() {
                      return $(this).attr('data-title');
                    },
                    container: 'body',
                    placement: 'top'
                  });
                };
                activateFgseaTooltips();
                table.on('draw.dt', activateFgseaTooltips);
                return table;
              "),
              options = list(
                columnDefs = list(
                  list(targets = c(0,1,2,3,4), visible = TRUE)
                )
              ),
              filter = 'top')
   
   pathwayTable <- DT::formatSignif(pathwayTable, columns = c("pval", "padj"), digits = 3)
   pathwayTable <- DT::formatRound(pathwayTable, columns = c("ES", "NES"), digits = 2)
   pathwayTable
   
  })

  output$pathwayAnalysisTopDotPlotHeading <- renderUI({
    tags$h3(paste0("Top 20 Pathways: ", input$selectIn2, " vs ", input$selectIn1))
  })

  output$pathwayAnalysisTopDotPlot <- renderPlot({
    pathwayPlotData <- as.data.frame(fgseaResults) %>%
      filter(!is.na(padj)) %>%
      mutate(
        ## constrain -log10 FDR so that outliers do not take up too much of the screen
        adjustedPValueLog10 = scales::squish(-log10(padj), range = c(1, 20))
      )

    shiny::validate(
      shiny::need(nrow(pathwayPlotData) > 0, "No Hallmark FGSEA pathways available.")
    )

    topUpregulated <- pathwayPlotData %>%
      filter(NES > 0) %>%
      slice_max(order_by = NES, n = 10, with_ties = FALSE) %>%
      mutate(.direction = "Upregulated")
    topDownregulated <- pathwayPlotData %>%
      filter(NES < 0) %>%
      slice_min(order_by = NES, n = 10, with_ties = FALSE) %>%
      mutate(.direction = "Downregulated")
    topPathwayPlotData <- bind_rows(topDownregulated, topUpregulated) %>%
      mutate(.direction = factor(.direction, levels = c("Downregulated", "Upregulated")))

    shiny::validate(
      shiny::need(nrow(topPathwayPlotData) > 0, "No upregulated or downregulated Hallmark FGSEA pathways available.")
    )

    nesLimit <- max(abs(topPathwayPlotData$NES), na.rm = TRUE)
    pathwayPlot <- ggplot(
      topPathwayPlotData,
      aes(
        x = NES,
        y = reorder(pathway, -NES),
        fill = NES,
        size = adjustedPValueLog10
      )
    ) +
      geom_segment(
        data = topPathwayPlotData,
        aes(
          x = -Inf,
          xend = NES,
          y = reorder(pathway, -NES),
          yend = reorder(pathway, -NES)
        ),
        color = "grey80",
        linewidth = 0.3,
        inherit.aes = FALSE
      ) +
      geom_point(shape = 21, color = "grey35", stroke = 0.3) +
      facet_wrap(~ .direction, ncol = 1, scales = "free") +
      scale_x_continuous(expand = expansion(mult = 0.05)) +
      scale_y_discrete(
        expand = expansion(add = 0.6),
        labels = function(pathway) {
          pathway %>%
            str_replace_all("_", " ") %>%
            str_trunc(width = 75)
        }
      ) +
      scale_fill_gradient2(
        name = "NES",
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        limits = c(-nesLimit, nesLimit)
      ) +
      scale_size_continuous(
        name = "-log10 adjusted p-value",
        range = c(3, 12),
        breaks = c(1, 2, 3),
        labels = c("1", "2", "3"),
        limits = c(1, 20)
      ) +
      labs(
        x = "Normalized Enrichment Score (suppressed <- 0 -> activated)",
        y = NULL
      ) +
      theme_classic() +
      theme(
        text = element_text(size = 12),
        axis.text.y = element_text(size = 10, margin = margin(r = 2)),
        strip.text = element_text(size = 13),
        panel.spacing.y = grid::unit(0.1, "lines")
      )

    pathwayPlot
  })
}

### DEG Calculation --------------------------------------------------------------------------------------------------
calculateWilcoxonDEG <- function(rnaSeqData1, rnaSeqData2) {
  validCellLines1 <- colSums(is.na(rnaSeqData1)) != nrow(rnaSeqData1)
  validCellLines2 <- colSums(is.na(rnaSeqData2)) != nrow(rnaSeqData2)
  
  rnaSeqData1 <- rnaSeqData1[, validCellLines1]
  rnaSeqData2 <- rnaSeqData2[, validCellLines2]

  # Perform Wilcoxon rank-sum test for each gene
  pValues <- sapply(1:nrow(rnaSeqData1), function(i) {
    wilcox.test(rnaSeqData1[i, ], rnaSeqData2[i, ], exact = FALSE)$p.value
  })

  # calculate FDR
  fdr <- p.adjust(pValues, method = "BH")

  # Average expression values across cell lines for each group
  avgExpr1 <- rowMeans(rnaSeqData1, na.rm = TRUE)
  avgExpr2 <- rowMeans(rnaSeqData2, na.rm = TRUE)

  # Calculate fold change
  foldChange <- avgExpr2 - avgExpr1

  # Create results table with fold change
  results <- data.frame(
    Gene = rownames(rnaSeqData1),
    Log2_Fold_Change = round(foldChange, 2),
    P_Value = round(pValues, 2),
    FDR = round(fdr, 2)
  )

  #Remove 'xsq' prefix from gene names
  results$Gene <- sub("^xsq", "", results$Gene)

  #Gene set annotation
  if (require(geneSetPathwayAnalysis)){
    results$Annotation <- geneSetPathwayAnalysis::geneAnnotTab[match(results$Gene,rownames(geneSetPathwayAnalysis::geneAnnotTab)), "SHORT_ANNOT"]
    results$Annotation[is.na(results$Annotation)] <- ""
  }

  results
}

calculateLimmaFit <- function(expressionData, group) {
  group <- factor(group, levels = unique(group))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- LETTERS[seq_len(nlevels(group))]

  groupValidValues <- lapply(levels(group), function(groupName) {
    groupData <- expressionData[, group == groupName, drop = FALSE]
    rowSums(is.finite(groupData)) > 0
  })

  ## Preprocessing: Filter out genes with NA or Inf values
  validValues <- Reduce(`&`, groupValidValues)

  ## Preprocessing: Filter out genes with 0 variance
  geneVar <- apply(expressionData, 1, var, na.rm = TRUE)
  variableGene <- !is.na(geneVar) & geneVar > 0

  ## Preprocessing: Apply Filters
  geneKeep <- validValues & variableGene
  expressionData <- expressionData[geneKeep, , drop = FALSE]
  
  fit <- limma::lmFit(expressionData, design)
  list(fit = fit, design = design, groupLevels = levels(group))
}

calculateLimmaContrast <- function(limmaFit, controlGroup, testGroup) {
  controlColumn <- colnames(limmaFit$design)[match(controlGroup, limmaFit$groupLevels)]
  testColumn <- colnames(limmaFit$design)[match(testGroup, limmaFit$groupLevels)]
  contrast <- limma::makeContrasts(
    contrasts = paste(testColumn, "-", controlColumn),
    levels = limmaFit$design
  )

  fit2 <- limma::contrasts.fit(limmaFit$fit, contrast)
  fit2 <- limma::eBayes(fit2, trend = TRUE)
  top <- limma::topTable(fit2, number = Inf, sort.by = "none")
  results <- data.frame(
    Gene = sub("^(xsq|exp)", "", rownames(top)),
    top,
    check.names = FALSE
  )
  
  if (require(geneSetPathwayAnalysis)){
    results$Annotation <- geneSetPathwayAnalysis::geneAnnotTab[match(results$Gene,rownames(geneSetPathwayAnalysis::geneAnnotTab)), "SHORT_ANNOT"]
    results$Annotation[is.na(results$Annotation)] <- ""
  }
  results
}
