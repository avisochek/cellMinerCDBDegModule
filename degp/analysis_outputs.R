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
