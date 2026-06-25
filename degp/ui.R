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
  tabsetPanel(
    tabPanel("Differential Expression Analysis",
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
  )
}
