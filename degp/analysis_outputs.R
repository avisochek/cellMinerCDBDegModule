renderAnalysisOutputs <- function(input, output, degResults, rnaSeqData1, rnaSeqData2, fgseaResults) {
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
              caption = HTML(paste("Differential Expression Analysis: ", input$selectIn1, " vs ", input$selectIn2,
                                   "<br><small>Search table by regular expressions</small>"))
              
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
    Significant=rep("Not Significant",length(log_FC))
    Significant[which(degResults$adj.P.Val<0.05 & degResults$logFC>=2)]="Upregulated: Log2 FoldChange > 2 & FDR < 0.05"
    Significant[which(degResults$adj.P.Val<0.05 & degResults$logFC<=-2)]="Downregulated: Log2 FoldChange < 2 & FDR < 0.05"
    
    gene = sub("^xsq", "", rownames(degResults))
    volcano_data=data.frame(gene,log_FC,log_pval,Significant)
    
    volcano_plot <- ggplot(volcano_data, aes(x = log_FC, y = log_pval, color = Significant, text = gene)) +
      geom_point() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
      scale_color_manual(values = c("blue", "grey", "red")) +
      labs(
        title = paste0('Volcano plot for: ', input$selectIn1, " vs ", input$selectIn2),
        x = "Log2 Fold Change",
        y = "-Log10 FDR"
      ) +
      scale_x_continuous(
        limits = c(-5, 5),
        breaks = seq(-5, 5, by = 0.5)
      ) +
      scale_y_continuous(
        limits = c(0, 10),
        breaks = seq(0, 10, by = 0.5)
      ) +
      theme_classic()
  
    volcano_plot <- ggplotly(volcano_plot, tooltip = "text")
    
    volcano_plot
    
  })
  
  
  ## Heatmap
  output$heatmapPlot <- renderPlotly({
    
    topUpregulated <- head(degResults[order(degResults$P.Value, -degResults$logFC), ], 10)
    topDownregulated <- head(degResults[order(degResults$P.Value, degResults$logFC), ], 10)
    
    topGenes <- rbind(topDownregulated, topUpregulated)
    geneNames <- c(rownames(topDownregulated), rownames(topUpregulated))
    #geneNames <- topGenes$Gene
    
    
    rnaSeqData1 <- data.frame(rnaSeqData1, check.names = FALSE) %>% select_if(~ !all(is.na(.)))
    rnaSeqData2 <- data.frame(rnaSeqData2, check.names = FALSE) %>% select_if(~ !all(is.na(.)))
    
    # print(dim(rnaSeqData1))
    # print(dim(rnaSeqData2))
    
    
    expressionData <- cbind(rnaSeqData1, rnaSeqData2)
    # print(dim(expressionData))
    
    # colnames(expressionData) <- c(paste(colnames(rnaSeqData1), "Group1", sep="_"),
    #                                       paste(colnames(rnaSeqData2), "Group2", sep="_"))
    
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
              scale = "row",
              # grid_gap = 1,
              fontsize_col = 5,
              fontsize_row = 5
              )
    
  })
  
  
  # Display FGSEA results
  output$pathwayAnalysisResults <- renderDT({
   pathwayTable <- datatable(fgseaResults,
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

  output$pathwayAnalysisDotPlot <- renderPlot({
    pathwayPlotData <- as.data.frame(fgseaResults) %>%
      filter(startsWith(pathway, "HALLMARK_"), !is.na(padj)) %>%
      mutate(
        .sign = if_else(NES >= 0, "activated", "suppressed"),
        leadingEdgeCount = lengths(leadingEdge)
      )

    shiny::validate(
      shiny::need(nrow(pathwayPlotData) > 0, "No Hallmark FGSEA pathways available.")
    )

    nesLimit <- max(abs(pathwayPlotData$NES), na.rm = TRUE)

    ggplot(pathwayPlotData, aes(x = NES, y = reorder(pathway, NES), color = .sign, size = leadingEdgeCount)) +
      geom_point() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
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

  output$pathwayAnalysisTopDotPlot <- renderPlot({
    pathwayPlotData <- as.data.frame(fgseaResults) %>%
      filter(startsWith(pathway, "HALLMARK_"), !is.na(padj)) %>%
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
