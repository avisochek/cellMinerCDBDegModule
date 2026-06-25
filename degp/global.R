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

source("degp/calculateDEG.R", local = TRUE)
source("degp/combineGeneSets.R", local = TRUE)
# Define gene sets for pathway analysis 
hallmarkGeneSets <- gmtPathways("degp/h.all.v2026.1.Hs.symbols.gmt")
dtbGeneSets <- gmtPathways("degp/220411.Genelist.gmt")
reactomeGeneSets <- gmtPathways("degp/c2.cp.reactome.v2026.1.Hs.symbols.gmt")

# Combine the gene sets
combinedGeneSets <- combineGeneSets(hallmarkGeneSets, dtbGeneSets)
combinedGeneSets <- combineGeneSets(combinedGeneSets, reactomeGeneSets)
