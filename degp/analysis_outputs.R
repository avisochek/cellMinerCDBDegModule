renderAnalysisOutputs <- function(input, output, degResults, rnaSeqData1, rnaSeqData2, fgseaResults) {
  # Render in results table
  output$resultsTable <- DT::renderDT({
    datatable(degResults, 
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
  })
  
  
  output$volcanoPlot <- renderPlotly({
    # # Generate volcano plot
    # degResults$Significance <- ifelse(degResults$P_Value < 0.05 & abs(degResults$Log2_Fold_Change) > 2, "Significant", "Not Significant")
    # 
    # y_limit <- max(-log10(degResults$P_Value)) * 1.5 # Extend limit by 50% for aesthetics
    # 
    # plot <- ggplot(degResults, aes(x = Log2_Fold_Change, y = -log10(P_Value), color = Significance)) +
    #   geom_point(size = 0.5, alpha = 0.5) +
    #   theme_minimal() +
    #   scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
    #   labs(title = "Volcano Plot of Differentially Expressed Genes",
    #        x = "Log2 Fold Change",
    #        y = "-Log10 P-value")
    # 
    # ggplotly(plot)
    
    degResults <- degResults[abs(degResults$Log2_Fold_Change) >= 0.10, ]
    degResults$P_Value[degResults$P_Value == 0] <- 10e-10
    
    
    log_FC = degResults$Log2_Fold_Change
    log_pval = round(-log10(degResults$P_Value), 2)
    Significant=rep("Not Significant",length(log_FC))
    Significant[which(degResults$P_Value<0.05 & abs(degResults$Log2_Fold_Change)>=2)]="Log2 FoldChange > +/- 2 & P-Value < 0.05"
    Significant[which(degResults$P_Value<0.05 & abs(degResults$Log2_Fold_Change)<1)]="P-Value < 0.05"
    Significant[which(degResults$P_Value>=0.05 & abs(degResults$Log2_Fold_Change)>=1)]="Log2 FoldChange > +/- 1"
    
    gene = degResults$Gene
    volcano_data=as.data.frame(cbind(gene,log_FC,log_pval,Significant))
    
    # Create the plotly object
    volcano_plot <- plot_ly(
      type = 'scatter',
      data = volcano_data,
      x = ~log_FC,
      y = ~log_pval,
      text = ~gene,
      #mode = "markers",
      color = ~Significant,
      colors = c("blue", "black", "red")
    )
    
    
    # Configure the layout of the plot
    volcano_plot <- layout(
      volcano_plot,
      title = paste0('Volcano plot for: ', input$selectIn1, " vs ", input$selectIn2),
      xaxis = list(
        title = "Log2 Fold Change",
        range = c(-5, 5),
        tickvals = seq(-5, 5, by = 0.5)
      ),
      yaxis = list(
        title = "-Log10 P-Value",
        range = c(0, 10),
        tickvals = seq(0, 10, by = 0.5)
      )
    )
    
    volcano_plot
    
  })
  
  
  ## Heatmap
  output$heatmapPlot <- renderPlotly({
    
    topUpregulated <- head(degResults[order(degResults$P_Value, -degResults$Log2_Fold_Change), ], 10)
    topDownregulated <- head(degResults[order(degResults$P_Value, degResults$Log2_Fold_Change), ], 10)
    
    topGenes <- rbind(topDownregulated, topUpregulated)
    geneNames <- paste0("xsq", topGenes$Gene)
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
   
   formatRound(pathwayTable, columns = c("pval", "padj", "ES", "NES"), digits = 2)
   
  })
}
