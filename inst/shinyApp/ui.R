shinyUI(fluidPage(
    
    titlePanel("IntLim (Integrate Metabolomics with other Omics data using linear modeling"),
    sidebarLayout(
        sidebarPanel(
            
            
                
            fileInput('file1', 'Choose CSV File',
                      accept=c('text/csv', 
                               'text/comma-separated-values,text/plain', 
                               '.csv')),
            textInput("metabid", "Metab ID", "BIOCHEMICAL"),
            textInput("geneid", "Gene ID", "X"),
            
            checkboxGroupInput(inputId="groupA",
                                   label="Pick Types for Group A", 
                                   choices=c("Normal"="N","Tumor"="T"), selected = "Tumor"),
                
           checkboxGroupInput(inputId="groupB",
                                   label="Pick Types for Group B", 
                                   choices=c("Normal"="N","Tumor"="T"), selected = "Tumor")
            
        ), # end sidebarPanel
            
        mainPanel(
            tabsetPanel(
                tabPanel('Data summary',
                       pre(textOutput('contents'))
                )
             ) # end tabsetPanel
               # tabPanel('Analyze By Group',pre(textOutput('group')))
                
            ) # end mainPanel
        ) # end sidebarLayout
    ) # end fluidPage
) # end shinUI

