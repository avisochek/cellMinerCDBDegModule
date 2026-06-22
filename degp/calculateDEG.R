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
  validCellLines <- colSums(is.na(expressionData)) != nrow(expressionData)

  expressionData <- expressionData[, validCellLines]
  group <- factor(group[validCellLines], levels = unique(group[validCellLines]))
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
