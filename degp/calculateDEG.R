#' PURPOSE_FIXME
#' 
#' @param rnaSeqData1 RNA-seq expression values from the control cell line
#' @param rnaSeqData2 RNA-seq expression values from the test cell line
#' 
#' @@examples
#' # example code matrix(runif(X))
#' 
#' @author Nana A. Kusi
#' 
calculateDEG <- function(rnaSeqData1, rnaSeqData2) {
  # #Create vectors to store DEG p-values and W-scores
  # pValues <- numeric(nrow(rnaSeqData1))
  # wScores <- numeric(nrow(rnaSeqData1))
  
  # Remove NAs
  
  validCellLines1 <- colSums(is.na(rnaSeqData1)) != nrow(rnaSeqData1)
  validCellLines2 <- colSums(is.na(rnaSeqData2)) != nrow(rnaSeqData2)
  
  rnaSeqData1 <- rnaSeqData1[, validCellLines1]
  rnaSeqData2 <- rnaSeqData2[, validCellLines2]
  
  # Perform Wilcoxon rank-sum test for each gene
  pValues <- sapply(1:nrow(rnaSeqData1), function(i) {
    wilcox.test(rnaSeqData1[i, ], rnaSeqData2[i, ], exact = FALSE)$p.value
  })
  
  # for(i in 1:nrow(rnaSeqData1)) {
  #   testResult <- wilcox.test(rnaSeqData1[i, ], rnaSeqData2[i, ], exact = FALSE)
  #   pValues[i] <- testResult$p.value
  #   wScores[i] <- testResult$statistic
  # }
  
  # calculate FDR
  fdr <- p.adjust(pValues, method = "BH")
  
  # Average expression values across cell lines for each group 
  avgExpr1 <- rowMeans(rnaSeqData1, na.rm = TRUE)
  avgExpr2 <- rowMeans(rnaSeqData2, na.rm = TRUE)
  
  print(head(avgExpr1))
  print(head(avgExpr2))
  
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