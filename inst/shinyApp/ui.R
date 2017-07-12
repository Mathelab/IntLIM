
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
            "adjusted p value",
            tabName = "adPval",
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
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "about",
                                shiny::tabPanel("About",
                                                box(
                             width = 12,
                             #includeMarkdown("include.md")
                             tags$b("Here is some introduction of the InLim"),
                             h5("We developed a new Bioconductor R pachage and R shiny app  Utilizing our novel linear modeling approach, we are able to identify putative changes in gene-metabolite associations resulting from cancer-type (with respect to NCI-60 data) and tissue-type (with respect to breast cancer data).  The gene-metabolite pairs identified were enriched for pathways linked to cancer in both cases.  The IntLiM package and App is designed such that other researchers can better interpret their metabolomics data in light of transcriptomic data by using our new linear modeling approach.  While this is not the “be all end all” tool for interpreting data, this package and approach will greatly assist researchers in formulating novel hypothesis and proposing new studies especially with regards to the gene-metabolite pairs identified.  Integrating the results with pathway analysis tools will provide further insight.  The IntLim package and App will be available for download via Github and a sample data-set and vignette is provided for users (https://github.com/mingrui-liu/IntLim/).  With metabolomics being a rising field, it will be necessary to interpret metabolomics data in light of gene expression.  IntLim provides a novel approach for doing so.   
")
                         )
                )
        ),
        shinydashboard::tabItem(tabName = "loaddata",
                fluidRow(
                    
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
                        textInput("metabid", "Metab ID", "id"),
                        textInput("geneid", "Gene ID", "id"),
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
                        fileInput2('file1', 'Choose CSV File',labelIcon = "folder-open-o",
                                  accept=c('text/csv', 
                                           'text/comma-separated-values,text/plain', 
                                           '.csv'),progress = FALSE),
                        verbatimTextOutput('filename'),
                        actionButton("run", "Run"),
                        hr(),
                        tags$b("The statistic summary of the data"),
                        pre(tableOutput('stats')),
                        hr(),
                        tags$b("Verify the distribution of the input data."),
                        uiOutput("plot")
                        
                        
                        
                            )
                    
                    
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
                        tags$b("Please chose the filtered options"),
                        checkboxInput('filter', 'FILTER', FALSE),
                        numericInput("geneperc", "percentile cutoff for filtering genes:", 15, min = 0, max = 100),
                        numericInput("metabperc", "percentile cutoff for filtering metabolites:", 15, min = 0, max = 100),
                        actionButton("run2", "Run"),
                        hr(),
                        tags$b("The statistic summary of filtered data"),
                        pre(tableOutput('Fstats')),
                        
                        hr(),
                        tags$b("Verify the distribution of the filtered data."),
                        uiOutput('Fplot'),
                        hr()
                        
                        
                    )
                )
        ),
        
        shinydashboard::tabItem(tabName = "adPval",
                fluidRow(
                    box(
                        title = strong("Adjusted P value") ,
                        width = NULL,
                        solidHeader = TRUE,
                        h5("Run the linear models and plot distribution of p-values:"),
                        radioButtons("dataset", label = h3("Data set"),
                                     choices = list("metabolite" = "metabolite", "gene" = "gene"), 
                                     selected = "metabolite"),
                        hr(),
                        uiOutput('choosestype'),
                        actionButton("run3", "Plot the adjusted p value distribution"),
                        pre(plotOutput('Pdist')),
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
                        actionButton("run4", "Run heatmap"),
                        hr(),
                        highcharter::highchartOutput("heatmap"),
                        hr()
                        
                        )
                    )
                ),
        shinydashboard::tabItem(tabName = "scatterplot",
                                fluidRow(
                                    box(
                                        title = strong("Scatter plot the data") ,
                                        width = NULL,
                                        solidHeader = TRUE,
                                        h5("scatter plot
                                           
                                           "),
                                        uiOutput('choosesampletype'),
                                        hr(),
                                        tags$b("Top 6 pairs of difference of correlation "),
                                        pre(tableOutput('table')),
                                      
                                        hr(),
                                        tags$b("Please input the gene name and metab name you choose"),
                                        
                                        uiOutput('chooseMetabName'),
                                        uiOutput('chooseGeneName'),
        
                                        
                                        hr(),
                                        actionButton("run5", "Run scatter plot"),
                                        hr(),
                                        highcharter::highchartOutput("scatterPlot"),
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
