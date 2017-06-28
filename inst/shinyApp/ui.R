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
            "Correlation",
            tabName = "correlation",
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
            actionButton("buttonstop", strong("Click to Exit Shiny App")),
            icon = icon("sign-out")
        )
    )
)

body <- shinydashboard::dashboardBody(
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "about",
                shiny::tabPanel("About",
                         box(
                             width = 12,
                             #includeMarkdown("include.md")
                             tags$b("Here is some introduction of the InLim")
                         )
                )
        ),
        shinydashboard::tabItem(tabName = "loaddata",
                fluidRow(
                    HTML("<div class='col-sm-4' style='min-width:
                         1000px !important;'>"),
                    box(
                        title = strong("Input menu file"),
                        width = NULL,
                        solidHeader = TRUE,
                        tags$b("Please be sure that all files noted in the CSV file,
                               including the CSV file, are in the same folder."),
                        h5("This step does the following: "),
                        tags$ul(
                            tags$li("Loads a metadata spreadsheet with a CSV file extention."),
                            tags$li("The metadata associated with data files to be analyzed in IntLim is supplied
                                    as a CSV file with two columns and 6 rows: 
                                    type,filenames
                                    metabData,myfilename
                                    geneData,myfilename
                                    metabMetaData,myfilename (optional)
                                    geneMetaData,myfilename (optional)
                                    sampleMetaData,myfilename"),
                            tags$li(" Note also that the input data files should be in a specific format:
                                    metabData: rows are metabolites, columns are samples
                                    geneData: rows are genes, columns are samples
                                    metabMetaData: rows are metabolites, features are columns
                                    geneMetaData: rows are genes, features https://stackoverflow.com/documentationare columns
                                    sampleMetaData: rows are samples, features are columns
                                    In addition, the first column of the sampleMetaData file is assumed to be the sample id, 
                                    and those sample ids should match the columns of metabData and geneData (e.g. it is required
                                    that all sample ids in the metabData and geneData are also in the sampleMetaDatafile)."),
                            tags$li("Prints out the statistic summary of the data.")
                            ),
                        
                        hr(),
                        tags$b("Please input the MetabID and the GeneID for your data "),
                        textInput("metabid", "Metab ID", "BIOCHEMICAL"),
                        textInput("geneid", "Gene ID", "X"),
                        hr(),
                        fileInput('file1', 'Choose CSV File',
                                  accept=c('text/csv', 
                                           'text/comma-separated-values,text/plain', 
                                           '.csv')),
                        actionButton("run", "Run"),
                        hr(),
                        tags$b("The statistic summary of the data"),
                        pre(tableOutput('stats')),
                        hr(),
                        tags$b("Verify the distribution of the input data."),
                        highchartOutput2("plot")
                        
                        
                       
                            ),
                    
                    HTML("</div>")
                    )
        ),
        
        
        shinydashboard::tabItem(tabName = "Filterdata",
                fluidRow(
                    box(
                        title = strong("Filter Data (optional)") ,
                        width = NULL,
                        solidHeader = TRUE,
                        h5("This step help you filter the data"),
                          
                            hr(),
                            checkboxInput('filter', 'FILTER', FALSE),
                            numericInput("geneperc", "percentile cutoff for filtering genes:", 15, min = 0, max = 100),
                            numericInput("metabperc", "percentile cutoff for filtering metabolites:", 15, min = 0, max = 100),
                            actionButton("run2", "Run"),
                            hr(),
                            pre(tableOutput('Fstats')),
                                
                            hr(),
                            #showOutput("Fplot","Highcharts"),
                            hr()
                                
                            
                            )
                        )
        ),
    
        shinydashboard::tabItem(tabName = "correlation",
                fluidRow(
                    box(
                        title = strong("correlation") ,
                        width = NULL,
                        solidHeader = TRUE,
                        h5("Run the linear models and plot distribution of p-values:"),
                        radioButtons("dataset", label = h3("Data set"),
                                     choices = list("metabolite" = "metabolite", "gene" = "gene"), 
                                     selected = "metabolite"),
                        hr(),
                        uiOutput('choosestype'),
                        pre(plotOutput("Pdist")),
                        hr()
                        
                        )
                        )
                    ),
        shinydashboard::tabItem(tabName = "processresult",
                fluidRow(
                    box(
                        title = strong("Process the result") ,
                        width = NULL,
                        solidHeader = TRUE,
                        h5("Process the results and filter pairs of genes-metabolites based on 
                            adjusted p-values and differences in correlation coefficients between groups 1 and 2.
                            Then plot heatmap of significant gene-metabolite pairs

                           "),
                        hr(),
                        highchartOutput("heatmap"),
                        hr()
                        
                    )
                )
        )
        
        
        
        )
                )






shinyUI(fluidPage(
    dashboardPage(
        headerbar,
        sidebar,
        body
    )
)
) # end shinUI
