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
            icon = icon("bullseye"),
            badgeLabel = "step 3"
        ),
        shinydashboard::menuItem(
            "Process result",
            tabName = "processresult",
            icon = icon("bullseye"),
            badgeLabel = "step 4"
        ),
        shinydashboard::menuItem(
            "Scatter plot",
            tabName = "scatterplot",
            icon = icon("bullseye"),
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
                        width = 8,
                        solidHeader = TRUE,
                        h5("This step takes a CSV file as input (see About) and performs the following:"),
                        tags$ul(
                             tags$li("Loads the CSV input file and checks that all files exist"),
                             tags$li("Reads in all files from the CSV input file and creates a MultiOmics object"),
                             tags$li("Outputs a statistic summary of the data loaded in")
                        )
                        ),
                    shinydashboard::box(
                        width = 4,
                        infoBoxOutput("statusbox1", width = NULL)
                    )
                ),#end of info flow
                fluidRow(
                    shinydashboard::box(
                        shinyFilesButton('file1',
                                         'Select CSV File',
                                         'Provide CSV File to Load Data',
                                         FALSE),
                        verbatimTextOutput('filename'),
                        uiOutput('idChooseM'),
                        uiOutput('idChooseG'),
                        hr(),
                        actionButton("run", "Run"),
                        
                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                         tags$div("Loading...",id="loadmessage"))
                    ),
                    shinydashboard::box(
                        tags$b("Summary Statistics"),
                        #plot.new(),
                        
                        pre(dataTableOutput('stats')),
                        tags$style(type="text/css", '#stats tfoot {display:none;}')
                    )
                ),#end of select file and stats flow
                fluidRow(
                    shinydashboard::box(
                        width = 12,
                        tags$b("Distribution of Input Data"),
                        pre(htmlOutput("plot"))
                    )
                       
                )#end of plot flow
                
        ), # end tab loaddata
        
        
        shinydashboard::tabItem(tabName = "Filterdata",
                fluidRow(
                    shinydashboard::box(
                        title = strong("Filter Data (optional)") ,
                        width = 8,
                        solidHeader = TRUE,
                        h5("This step allows you to filter the metabolomics or gene expression data by a user-defined percentile cutoff."),
                        hr(),
                        numericInput("geneperc", "percentile cutoff for filtering genes (0-1):", 0, min = 0, max = 1),
                        numericInput("metabperc", "percentile cutoff for filtering metabolites (0-1):", 0, min = 0, max = 1),
                        numericInput("metabmiss", "missing value percent cutoff for filtering metabolites (0-1)", 0,min=0,max=1),
                        actionButton("run2", "Run")
                    ),
                    shinydashboard::box(
                        width=4,
                        infoBoxOutput("statusbox2", width = NULL)
                    )
                    ), # end filter option flow
                fluidRow(
                    shinydashboard::box(
                        tags$b("The statistic summary of origin data"),
                        dataTableOutput('Ostats'),
                        tags$style(type="text/css", '#Ostats tfoot {display:none;}')
                    ),
                    shinydashboard::box(
                        tags$b("The statistic summary of filtered data"),
                        dataTableOutput('Fstats'),
                        tags$style(type="text/css", '#Fstats tfoot {display:none;}')
                    )
                ),#end stats comparison flow    
                fluidRow(
                    shinydashboard::box(
                        tags$b("The distribution of the origin data."),
                        uiOutput('Oplot')
                    ),
                    shinydashboard::box(
                        tags$b("Verify the distribution of the filtered data."),
                        uiOutput('Fplot')
                    )
                )#end plot comparison flow
                
        ),
        
        shinydashboard::tabItem(tabName = "RunLM",
                fluidRow(
                    
                    shinydashboard::box(
                        title = strong("Run IntLIM") ,
                        width = 8,
                        solidHeader = TRUE,
                        h5("This step performs the linear models for all combinations of gene:metabolite pairs and then plots distribution of p-values.  "),
                        radioButtons("dataset", label = h5("Select the outcome set:"),
                                     choices = list("metabolite" = "metabolite", "gene" = "gene"), 
                                     selected = "metabolite"),
                        hr(),
                        uiOutput('choosestype'),
                        actionButton("run3", "Run")
                        
                    ),
                    shinydashboard::box(
                        width = 4,
                        infoBoxOutput("statusbox3", width = NULL)
                    )
                ),
                fluidRow(
                    shinydashboard::box(
                        width = NULL,
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
                                     tags$div("Loading...(It can take several minutes, please be patient)",id="loadmessage")),
                    tags$b("Histogram of p-values."),
                    highcharter::highchartOutput("Pdist")
                )
                )
        ),
        shinydashboard::tabItem(tabName = "processresult",
                fluidRow(
                    shinydashboard::box(
                        title = strong("Process the result") ,
                        width = 8,
                        solidHeader = TRUE,
                        h5("Process the results and filter pairs of genes-metabolites based on 
                           adjusted p-values and differences in correlation coefficients between groups 1 and 2.
                           Then plot heatmap of significant gene-metabolite pairs
                           
                           "),
                        actionButton("run4", "Run")
                        ),
                    shinydashboard::box(
                        width = 4,
                        infoBoxOutput("statusbox4", width = NULL)
                    )
                    ),#end of info flow
                fluidRow(
                    shinydashboard::box(
                        width = NULL,
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
                    plotly::plotlyOutput("heatmap")
                    )
                )
                ),
        shinydashboard::tabItem(tabName = "scatterplot",
                                fluidRow(
                                    shinydashboard::box(
                                        title = strong("Scatter plot") ,
                                        width = 8,
                                        solidHeader = TRUE,
                                        h5("This step present the table of gene-metabolite pairs and the absolute value of their 
                                           correlation differece"),
                                        h5("You can plot the scatter plot of prefered gene-metabolite pairs by clicking table")
                                        ),
                                    shinydashboard::box(
                                        width = 4,
                                        infoBoxOutput("statusbox5", width = NULL)
                                    )
                                ),# end of info flow
                                fluidRow(
                                    shinydashboard::box(
                                        width = NULL,
                                        tags$b("Pairs of difference of correlation "),
                                        pre(DT::dataTableOutput('table')),
                                        
                                        hr(),
                                        actionButton("run5", "Run"),
                                        
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
                                        highcharter::highchartOutput("scatterPlot")
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
