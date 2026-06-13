### User Interface --------------------------------------------------------------------------------------------------

ui <- fluidPage(
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
                   paste("<label class='dataset' for='dataSet'>Select Dataset</label>",
                         "<select id='dataSet'>", options, "</select>")
                 ),
                 br(),
                 br(),
                 h4("Please select cell lines by clicking on a row in the table"),
                 br(),
                 h4("Or search by tissue type"),
                 uiOutput("tissueSelector"),
                 textInput("groupName", "Group Name", value = ""),
                 actionButton(inputId = "createGroup", label = "Add to Group"),
                 br(),
                 br(),
                 uiOutput("groupInfoDisplay"),
                 br(),
                 uiOutput("choice1"),
                 uiOutput("choice2"),
                 actionButton("runAnalysis", "Run Analysis"),
                 #uiOutput("groupNamesDisplay"), 
                 actionButton("resetSelections", "Reset Selections")
               ),
               mainPanel(
                 tabsetPanel(id = "mainTabset",
                             tabPanel("Input",
                                      fluidRow(
                                        DT::DTOutput("dataSetTable")
                                      )
                             ),
                             tabPanel("Results",
                                      fluidRow(
                                        withSpinner(DT::DTOutput("resultsTable")), 
                                        downloadButton("downloadResults", "Download Results")
                                      )),
                             tabPanel("Volcano Plot",
                                      plotlyOutput("volcanoPlot")), 
                             tabPanel("Heatmap",
                                      plotlyOutput("heatmapPlot")),
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
                                      DT::DTOutput("pathwayAnalysisResults"))
                 )
               )
             )
    )
  )
)
