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
                        numericInput("geneperc", "percentile cutoff (0-1) for filtering genes (e.g. remove genes with mean values < cutoff):", 0, min = 0, max = 1),
                        numericInput("metabperc", "percentile cutoff (0-1) for filtering metabolites (e.g. remove metabolites with mean values < cutoff):", 0, min = 0, max = 1),
                        numericInput("metabmiss", "missing value percent cutoff (0-1) for filtering metabolites (e.g. metabolites with > % cutoff missing values will be removed)", 0,min=0,max=1),
                        actionButton("run2", "Run")
                    ),
                    shinydashboard::box(
                        width=4,
                        infoBoxOutput("statusbox2", width = NULL),
                        downloadButton('downloadFdata', 'Download')
                        
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
                        width = 6,
                        height =650,
                        solidHeader = TRUE,
                        h5("This step performs the linear models for all combinations of gene:metabolite pairs and then plots distribution of p-values."),
			h5("The linear model performed is 'm ~ g + p + g:p' where "),
			tags$ul(
				tags$li("'m' is the metabolite abundance"),
				tags$li("‘g’ is the gene expression level"), 
				tags$li("‘p’ is the phenotype (e.g. tumor vs non-tumor)"),
				tags$li("‘g:p’ is the interaction between phenotype and gene expression")
			),
			h5("A statistically significant p-value of the the interaction term ‘g:p’ indicates that the gene-metabolite relationship is phenotype-specific. Please see manuscript for more details."),
#                        radioButtons("dataset", label = h5("Select column that has the categories that you wish to compare (e.g. tumor vs. non-tumor):"),
#                                     choices = list("metabolite" = "metabolite", "gene" = "gene"), 
#                                     selected = "metabolite"),
                        hr(),
                        uiOutput('choosestype'),
                        numericInput("nrpoints", "number of points to be plotted in lowest density areas:", 10000, min = 0, max = 30000),
                        numericInput("pvalcutoff1","cutoff of FDR-adjusted p-value for filtering(0 - 1) :", 0.05, min = 0, max = 1),
                        numericInput("diffcorr1", "cutoff of differences in correlations for filtering (0-1):", 0.5, min = 0, max = 1),
                        actionButton("run3", "Run")

                       
                        
                    ),
                    shinydashboard::box(
                        width = 6,
                        height =650,
                        
                       
                        
                        numericInput("breaks","Breaks of histogram",100,min=10,max = 500),
                        plotOutput("Pdist"),
                        hr(),
                        textOutput("Ptext")
                        
                    )
                ),
                fluidRow(
                    shinydashboard::box(
                        width = 8,
                        
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
                                     tags$div("Loading...(It might take several minutes depending on the size of the dataset,please be patient!)",id="loadmessage")),
                    
                    plotOutput("volcanoPlot")
                    
                    
                ),
                shinydashboard::box(
                    width = 4,
                    infoBoxOutput("statusbox3", width = NULL)
                    
                
                )
                )
        ),
        shinydashboard::tabItem(tabName = "processresult",
                fluidRow(
                    shinydashboard::box(
                        title = strong("Process the result") ,
                        width = 8,
                        height = 200,
                        solidHeader = TRUE,
                        h5("Process the results and filter pairs of genes-metabolites based on 
                           adjusted p-values and differences in correlation coefficients between groups 1 and 2.
                           Then plot heatmap of significant gene-metabolite pairs
                           
                           ")
                        
                        ),
                    shinydashboard::box(
                        width = 4,
                        height = 200,
                        infoBoxOutput("statusbox4", width = NULL),
                        downloadButton('downloadData', 'Download')
                    )
                    ),#end of info floww
                fluidRow(
                    shinydashboard::box(
                        width = 4,
                        numericInput("pvalcutoff","cutoff of FDR-adjusted p-value for filtering(0 - 1) :",NA, min = 0, max = 1),
                        numericInput("diffcorr", "cutoff of differences in correlations for filtering (0-1):",NA, min = 0, max = 1),
                        textInput("corrtype","spearman or pearson or other parameters allowed by cor() function","spearman"),
                        actionButton("run4", "Run")
                    ),
                    
                    shinydashboard::box(
                        width = 8,
                        
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
                                        
                                        #pre(textOutput("temp")),
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
                                                         tags$div("Loading...",id="loadmessage"))
                                        
                                        
                                    )
                                ),#end of table flow
                                    fluidRow(
                                    shinydashboard::box(
                                        width = NULL,
                                        pre(uiOutput("scatterplot"))
                                    
                                    )
                                    )#end of scatterplot flow
        
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
