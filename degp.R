# Load packages ----------------------------------------------------------------
library(shiny)
library(rcellminer)
library(shinycssloaders)
library(plotly)
library(DT)
library(jsonlite)
library(stringr)
library(fgsea)
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

# Define gene sets for pathway analysis 
hallmarkGeneSets <- gmtPathways("h.all.v2026.1.Hs.symbols.gmt")
dtbGeneSets <- gmtPathways("220411.Genelist.gmt")
reactomeGeneSets <- gmtPathways("c2.cp.reactome.v2026.1.Hs.symbols.gmt")

# Combine the gene sets
combinedGeneSets <- combineGeneSets(hallmarkGeneSets, dtbGeneSets)
combinedGeneSets <- combineGeneSets(combinedGeneSets, reactomeGeneSets)

### User Interface --------------------------------------------------------------------------------------------------
degpInput <- function(id) {
  ns <- NS(id)
  tabPanel("Differential Expression Analysis",
  fluidPage(
  useShinyjs(),
    tags$head(
      tags$style(HTML("
        .shiny-notification {
          position: fixed !important;
          top: 40% !important;
          left: 50% !important;
          width: auto !important;
          padding: 10px !important;
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
                 uiOutput(ns("tissueSelector")),
                 actionButton(inputId = ns("addGroup"), label = "Add Group"),
                 br(),
                 br(),
                 uiOutput(ns("groupInfoDisplay")),
                 br(),
                 tags$div(class = "sidebar-step-heading", "3. Select Contrast"),
                 selectizeInput(ns("selectIn1"), "Select Control Group:", choices = c("")),
                 selectizeInput(ns("selectIn2"), "Select Test Group:", choices = c("")),
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
                                      fluidRow(
                                        withSpinner(DT::DTOutput(ns("resultsTable"))),
                                        downloadButton(ns("downloadResults"), "Download Results")
                                      )),
                             tabPanel("Volcano Plot",
                                      plotlyOutput(ns("volcanoPlot"))),
                             tabPanel("Heatmap",
                                      plotlyOutput(ns("heatmapPlot"))),
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
                                                 DT::DTOutput(ns("pathwayAnalysisResults"))),
                                        tabPanel("Top FGSEA Plot",
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
    dataSourceChoices <- setNames(
      names(srcContent),
      vapply(expressionConfig[names(srcContent)], function(x) { x[["displayName"]] }, character(1))
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
    data <- req(state$sampleData())
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
  
  
  
}
### Analysis Outputs --------------------------------------------------------------------------------------------------
renderAnalysisOutputs <- function(input, output, degResults, exprData1, exprData2, fgseaResults) {
  # Render in results table
  output$resultsTable <- DT::renderDT({
    resultsTable <- datatable(degResults, 
              options = list(
                search = list(regex = TRUE),
                columnDefs = list(
                  list(visible = FALSE, targets = 0) 
                )
              ),
              filter = 'top',
              caption = HTML(paste0(
                "<div style='font-size: 24px; line-height: 1.35;'>",
                "Differential Expression Analysis: ", input$selectIn1, " vs ", input$selectIn2,
                "<br><small style='font-size: 13px;'>Search table by regular expressions</small>",
                "</div>"
              ))
              
    )
    
    resultsTable <- DT::formatRound(resultsTable, columns = c("logFC", "AveExpr", "t", "B"), digits = 2)
    resultsTable <- DT::formatSignif(resultsTable, columns = c("P.Value", "adj.P.Val"), digits = 3)
    resultsTable
  })
  
  
  output$volcanoPlot <- renderPlotly({
    
    degResults <- degResults[abs(degResults$logFC) >= 0.10, ]
    degResults$P.Value[degResults$adj.P.Val == 0] <- 10e-10
    
    
    log_FC = degResults$logFC
    log_pval = -log10(degResults$adj.P.Val)
    upregulatedLabel <- "Upregulated: Log2 FoldChange >= 1 & FDR < 0.05"
    downregulatedLabel <- "Downregulated: Log2 FoldChange <= -1 & FDR < 0.05"
    topGeneLabel <- "Top 10 upregulated / downregulated"
    
    Significant=rep("Not Significant",length(log_FC))
    Significant[which(degResults$adj.P.Val<0.05 & degResults$logFC>=1)]=upregulatedLabel
    Significant[which(degResults$adj.P.Val<0.05 & degResults$logFC<=-1)]=downregulatedLabel
    
    gene = sub("^xsq", "", rownames(degResults))
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
    topGeneData <- volcano_data[!is.na(volcano_data$TopGene), ]
    pointSize <- 1.5
    xAxisMax <- max(2, abs(log_FC), na.rm = TRUE)
    yAxisMax <- max(2, log_pval, na.rm = TRUE)
    
    volcano_plot <- ggplot(volcano_data, aes(x = log_FC, y = log_pval, color = Significant, text = gene)) +
      geom_point(size = pointSize) +
      geom_point(
        data = topGeneData,
        aes(color = TopGene),
        shape = 1,
        fill = NA,
        size = pointSize,
        stroke = 0.5
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
      scale_color_manual(values = c(
        "Upregulated: Log2 FoldChange >= 1 & FDR < 0.05" = "red",
        "Downregulated: Log2 FoldChange <= -1 & FDR < 0.05" = "blue",
        "Not Significant" = "grey",
        "Top 10 upregulated / downregulated" = "gold"
      )) +
      labs(
        title = paste0('Volcano plot for: ', input$selectIn1, " vs ", input$selectIn2),
        x = "Log2 Fold Change",
        y = "-Log10 FDR"
      ) +
      scale_x_continuous(
        limits = c(-xAxisMax, xAxisMax)
      ) +
      scale_y_continuous(
        limits = c(0, yAxisMax)
      ) +
      theme_classic()
  
    volcano_plot <- ggplotly(volcano_plot, tooltip = "text")
    
    volcano_plot
    
  })
  
  
  ## Heatmap
  output$heatmapPlot <- renderPlotly({
    
    topUpregulated <- head(degResults[order(-degResults$logFC, degResults$P.Value), ], 10)
    topDownregulated <- head(degResults[order(degResults$logFC, degResults$P.Value), ], 10)
    
    geneNames <- c(rownames(topDownregulated), rownames(topUpregulated))
    #geneNames <- topGenes$Gene
    
    
    exprData1 <- data.frame(exprData1, check.names = FALSE) %>% select_if(~ !all(is.na(.)))
    exprData2 <- data.frame(exprData2, check.names = FALSE) %>% select_if(~ !all(is.na(.)))
    
    # print(dim(exprData1))
    # print(dim(exprData2))
    
    
    expressionData <- cbind(exprData1, exprData2)
    # print(dim(expressionData))
    
    # colnames(expressionData) <- c(paste(colnames(exprData1), "Group1", sep="_"),
    #                                       paste(colnames(exprData2), "Group2", sep="_"))
    
    expressionData <- expressionData[geneNames, , drop = FALSE]
    #Remove 'xsq' prefix from gene names
    rownames(expressionData) <- gsub("^xsq", "", rownames(expressionData))
    
    heatmaply(expressionData, 
              main = paste0("Heatmap of top 10 upregulated genes and bottom 10 downregulated genes in ", input$selectIn1, " vs ", input$selectIn2),
              xlab = "Cell Lines", 
              ylab = "Genes", 
              colors = colorRampPalette(c("navy", "white", "firebrick"))(255),
              show_rownames = TRUE,
              show_colnames = TRUE,
              show_dendrogram = c(FALSE, FALSE),
              margins = c(NA, NA, 150, NA), ## Needed to make room for title
              scale = "row",
              # grid_gap = 1,
              fontsize_col = 5,
              fontsize_row = 5
              )
    
  })
  
  
  # Display FGSEA results
  output$pathwayAnalysisResults <- renderDT({
   truncateFgseaCell <- function(value) {
     value <- paste(value, collapse = ", ")
     label <- htmltools::htmlEscape(value)
     tooltip <- htmltools::htmlEscape(gsub("_", "_ ", value), attribute = TRUE)
     sprintf("<div style='max-width: 260px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;' data-toggle='tooltip' data-title='%s'>%s</div>", tooltip, label)
   }
   
   fgseaResults$pathway <- vapply(fgseaResults$pathway, truncateFgseaCell, character(1))
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
              filter = 'top',
              caption = paste("Pathway Analysis Results: ", input$selectIn1, " vs ", input$selectIn2))
   
   pathwayTable <- DT::formatSignif(pathwayTable, columns = c("pval", "padj"), digits = 3)
   pathwayTable <- DT::formatRound(pathwayTable, columns = c("ES", "NES"), digits = 2)
   pathwayTable
   
  })

  output$pathwayAnalysisTopDotPlot <- renderPlot({
    pathwayPlotData <- as.data.frame(fgseaResults) %>%
      filter(!is.na(padj)) %>%
      mutate(
        .sign = if_else(NES >= 0, "activated", "suppressed"),
        leadingEdgeCount = lengths(leadingEdge)
      )

    shiny::validate(
      shiny::need(nrow(pathwayPlotData) > 0, "No Hallmark FGSEA pathways available.")
    )

    nesLimit <- max(abs(pathwayPlotData$NES), na.rm = TRUE)
    topUpregulated <- pathwayPlotData %>%
      filter(NES > 0) %>%
      slice_max(order_by = NES, n = 10, with_ties = FALSE) %>%
      mutate(.direction = "Top 10 upregulated")
    topDownregulated <- pathwayPlotData %>%
      filter(NES < 0) %>%
      slice_min(order_by = NES, n = 10, with_ties = FALSE) %>%
      mutate(.direction = "Top 10 downregulated")
    topPathwayPlotData <- bind_rows(topDownregulated, topUpregulated)

    shiny::validate(
      shiny::need(nrow(topPathwayPlotData) > 0, "No upregulated or downregulated Hallmark FGSEA pathways available.")
    )

    ggplot(topPathwayPlotData, aes(x = NES, y = reorder(pathway, NES), color = .sign, size = leadingEdgeCount)) +
      geom_point() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      facet_wrap(~ .direction, nrow = 1, scales = "free_y") +
      scale_x_continuous(limits = c(-nesLimit, nesLimit)) +
      scale_color_manual(name = "Direction", values = c("activated" = "red", "suppressed" = "blue")) +
      scale_size_continuous(name = "Leading edge") +
      labs(
        title = paste(input$selectIn1, "vs", input$selectIn2),
        x = "Normalized Enrichment Score (suppressed <- 0 -> activated)",
        y = NULL
      ) +
      theme_classic()
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
