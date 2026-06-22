### User Interface --------------------------------------------------------------------------------------------------

degpUI <- function(id) {
  ns <- NS(id)
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
    "))
  ),
  tabsetPanel(
    tabPanel("Differential Expression Analysis",
             sidebarLayout(
               #tags$head(tags$style(".shiny-notification {position: fixed; top: 60% ;left: 50%")),
               sidebarPanel(
                 HTML(
                   paste0("<label class='dataset' for='", ns("dataSet"), "'>Select Dataset</label>",
                          "<select id='", ns("dataSet"), "'>", options, "</select>")
                 ),
                 br(),
                 h2("Select Groups for DEG and Pathway Analysis:"),
                 h4("Either select cell lines from the table, or choose 1 or more tissue types from the 'Select Tissue' selector."),
                 h4("When you are done, click 'Add Group'."),
                 br(),
                 uiOutput(ns("tissueSelector")),
                 textInput(ns("groupName"), "Group Name", value = ""),
                 actionButton(inputId = ns("addGroup"), label = "Add Group"),
                 br(),
                 br(),
                 uiOutput(ns("groupInfoDisplay")),
                 br(),
                 selectizeInput(ns("selectIn1"), "Select Control Group:", choices = c("")),
                 selectizeInput(ns("selectIn2"), "Select Test Group:", choices = c("")),
                 actionButton(ns("runAnalysis"), "Run Analysis"),
                 actionButton(ns("reset"), "Reset")
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
  )
}
