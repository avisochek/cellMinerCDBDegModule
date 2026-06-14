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

### --------------------------------------------------------------------------------------------------
### LOAD CONFIGURATION AND REQUIRED DATA SOURCE PACKAGES.
### --------------------------------------------------------------------------------------------------
config <- jsonlite::fromJSON("config.json")

source("calculateDEG.R", local = TRUE)
source("combineGeneSets.R", local = TRUE)

dataSourceChoices <- setNames(names(config),
                              vapply(config, function(x) { x[["displayName"]] },
                                     character(1)))
options = "";

for(y in 1:length(dataSourceChoices)){

  if (dataSourceChoices[y]=="nci60")
  {
    options =  paste0(options,"<option value=",dataSourceChoices[y]," selected>",names(dataSourceChoices)[y],"</option>;")
  }
  else
  {
    options =  paste0(options,"<option value=",dataSourceChoices[y],">",names(dataSourceChoices)[y],"</option>;");
  }
}
srcContent = readRDS("srcContent.rds")

# Define gene sets for pathway analysis 
hallmarkGeneSets <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
dtbGeneSets <- gmtPathways("220411.Genelist.gmt")

# Combine the gene sets
combinedGeneSets <- combineGeneSets(hallmarkGeneSets, dtbGeneSets)
