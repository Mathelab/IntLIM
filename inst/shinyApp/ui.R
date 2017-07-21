headerbar <- shinydashboard::dashboardHeader(
    title = "IntLim",
    titleWidth = 270,
    shinydashboard::dropdownMenu(
        type = "notifications",
        shinydashboard::notificationItem(
            text = "Plots might take some time to display",
            icon("truck"),
            status = "warning"
        )
    )
)

sidebar <- shinydashboard::dashboardSidebar(
    width = 270,
    shinydashboard::sidebarMenu(
        shinydashboard::menuItem("About",
                 tabName = "about",
                 icon = icon("info")),
        shinydashboard::menuItem(
            "Load Data",
            tabName = "loaddata",
            icon = icon("folder-open"),
            badgeLabel = "step 1"
        ),
        shinydashboard::menuItem(
            "Filter Data (optional)",
            tabName = "Filterdata",
            icon = icon("bullseye"),
            badgeLabel = "step 2"
        ),
        shinydashboard::menuItem(
            "Run Linear Models",
            tabName = "RunLM",
            icon = icon("bolt"),
            badgeLabel = "step 3"
        ),
        shinydashboard::menuItem(
            "Process result",
            tabName = "processresult",
            icon = icon("pie-chart"),
            badgeLabel = "step 4"
        ),
        shinydashboard::menuItem(
            "Scatter plot",
            tabName = "scatterplot",
            icon = icon("star"),
            badgeLabel = "step 5"
        ),
        shinydashboard::menuItem(
            actionButton("buttonstop", strong("Click to Exit Shiny App")),
            icon = icon("sign-out")
        )
    )
)

body <- shinydashboard::dashboardBody(
    
    shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "about",
                             shiny::tabPanel("About",
                                             shinydashboard::box(
                             width = 12,
                             includeMarkdown("README.md")
                         )
                )
        ),
        shinydashboard::tabItem(tabName = "loaddata",
                fluidRow(
                    
                    shinydashboard::box(
                        title = strong("Load Data"),
                        width = NULL,
                        solidHeader = TRUE,
                        h5("This function takes a CSV file as input (see About) and performs the following:"),
                        tags$ul(
                             tags$li("Loads the CSV input file and checks that all files exist"),
                             tags$li("Reads in all files from the CSV input file and creates a MultiOmics object"),
                             tags$li("Outputs a statistic summary of the data loaded in")
                        ),
                        hr(),
                        tags$head(tags$style(HTML(
                            ".fileinput_2 {
                            width: 0.1px;
                            height: 0.1px;
                            opacity: 0;
                            overflow: hidden;
                            position: absolute;
                            z-index: -1;
                            }"
                       ))),
                        shinyFilesButton('file1',
		                 'Select CSV File',
		                 'Provide CSV File to Load Data',
                                  FALSE),
                        verbatimTextOutput('filename'),
		                uiOutput('idChooseM'),
		                uiOutput('idChooseG'),
                        hr(),
                        actionButton("run", "Run"),
                        hr(),
                        tags$head(tags$style(type="text/css", "
                              loadmessage {
                                             position: fixed;
                                             top: 0px;
                                             left: 0px;
                                             width: 100%;
                                             padding: 5px 0px 5px 0px;
                                             text-align: center;
                                             font-weight: bold;
                                             font-size: 100%;
                                             color: #000000;
                                             background-color: #CCFF66;
                                             z-index: 105;
                                             }
                                             ")),
                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                         tags$div("Loading...",id="loadmessage")),
                        tags$b("The statistic summary of the data"),
                        plot.new(),
		                
                        pre(dataTableOutput('stats')),
		                tags$style(type="text/css", '#stats tfoot {display:none;}'),
                        hr(),
                        tags$b("Verify the distribution of the input data."),
                        plot.new(),
                        uiOutput("plot")
                            )
                        )
                        ), # end tab loaddata
        
        
        shinydashboard::tabItem(tabName = "Filterdata",
                fluidRow(
                    shinydashboard::box(
                        title = strong("Filter Data (optional)") ,
                        width = NULL,
                        solidHeader = TRUE,
                        h5("This step allows you to filter the metabolomics or gene expression data by a user-defined percentile cutoff."),
                        hr(),
                        #checkboxInput('filter', 'FILTER', FALSE),
                        numericInput("geneperc", "percentile cutoff for filtering genes:", 0, min = 0, max = 100),
                        numericInput("metabperc", "percentile cutoff for filtering metabolites:", 0, min = 0, max = 100),
                        actionButton("run2", "Run"),
                        
                        pre(textOutput("temp2")),
                        hr(),
                        verbatimTextOutput('FiltMessage'),
                        
                        tags$b("The statistic summary of filtered data"),
                        pre(dataTableOutput('Fstats')),
                        tags$style(type="text/css", '#Fstats tfoot {display:none;}'),
                       
                        hr(),
                        tags$b("Verify the distribution of the filtered data."),
                        uiOutput('Fplot'),
                        hr()
                        
                        
                    )
                )
        ),
        
        shinydashboard::tabItem(tabName = "adPval",
                fluidRow(
                    
                    shinydashboard::box(
                        title = strong("Adjusted P value") ,
                        width = NULL,
                        solidHeader = TRUE,
                        h5("Run the linear models and plot distribution of p-values:"),
                        radioButtons("dataset", label = h3("Data set"),
                                     choices = list("metabolite" = "metabolite", "gene" = "gene"), 
                                     selected = "metabolite"),
                        hr(),
                        uiOutput('choosestype'),
                        actionButton("run3", "Run"),
                        
                        tags$head(tags$style(type="text/css", "
                        loadmessage {
                                             position: fixed;
                                             top: 0px;
                                             left: 0px;
                                             width: 100%;
                                             padding: 5px 0px 5px 0px;
                                             text-align: center;
                                             font-weight: bold;
                                             font-size: 100%;
                                             color: #000000;
                                             background-color: #CCFF66;
                                             z-index: 105;
                                             }
                                             ")),
                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                         tags$div("Loading...",id="loadmessage")),
                    
                        pre(plotOutput('Pdist')),
                        hr()
                        
                    )
                )
        ),
        shinydashboard::tabItem(tabName = "processresult",
                fluidRow(
                    shinydashboard::box(
                        title = strong("Process the result") ,
                        width = NULL,
                        solidHeader = TRUE,
                        h5("Process the results and filter pairs of genes-metabolites based on 
                           adjusted p-values and differences in correlation coefficients between groups 1 and 2.
                           Then plot heatmap of significant gene-metabolite pairs
                           
                           "),
                        actionButton("run4", "Run heatmap"),
                        hr(),
                        tags$head(tags$style(type="text/css", "
                        loadmessage {
                                             position: fixed;
                                             top: 0px;
                                             left: 0px;
                                             width: 100%;
                                             padding: 5px 0px 5px 0px;
                                             text-align: center;
                                             font-weight: bold;
                                             font-size: 100%;
                                             color: #000000;
                                             background-color: #CCFF66;
                                             z-index: 105;
                                             }
                                             ")),
                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                         tags$div("Loading...",id="loadmessage")),
                        highcharter::highchartOutput("heatmap"),
                        hr()
                        
                        )
                    )
                ),
        shinydashboard::tabItem(tabName = "scatterplot",
                                fluidRow(
                                    shinydashboard::box(
                                        title = strong("Scatter plot") ,
                                        width = NULL,
                                        solidHeader = TRUE,
                                        h5("some introduction of scatter plot
                                           
                                           "),
                                        uiOutput('choosesampletype'),
                                        hr(),
                                        tags$b("Top 6 pairs of difference of correlation "),
                                        pre(tableOutput('table')),
                                      
                                        hr(),
                                        tags$b("Input the gene name and metab name you choose if meta
                                           data files are provided for genes and/or metabolites."),
                                        hr(),
                                        uiOutput('chooseMetabName'),
                                        uiOutput('chooseGeneName'),
        
                                        
                                        hr(),
                                        actionButton("run5", "Plot"),
                                        
                                        tags$head(tags$style(type="text/css", "
                        loadmessage {
                                                             position: fixed;
                                                             top: 0px;
                                                             left: 0px;
                                                             width: 100%;
                                                             padding: 5px 0px 5px 0px;
                                                             text-align: center;
                                                             font-weight: bold;
                                                             font-size: 100%;
                                                             color: #37649b;
                                                             background-color: #CCFF66;
                                                             z-index: 105;
                                                             }
                                                             ")),
                                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                         tags$div("Loading...",id="loadmessage")),
                                        highcharter::highchartOutput("scatterPlot"),
                                        hr()
                                        
                                        )
                                    )
                                )
        
        
        
        
        
        )
                )





shinyUI(fluidPage(
    shinydashboard::dashboardPage(
        headerbar,
        sidebar,
        body
    )
)
) # end shinUI
