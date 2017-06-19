
library(shiny)
library(MultiAssayExperiment)
source('~/Documents/BIM/IntLim-master/IntLim/R/MetaboliteSet_addMetabolite.R')
source('~/Documents/BIM/IntLim-master/IntLim/R/internalfunctions.R')
source('~/Documents/BIM/IntLim-master/IntLim/R/AllClasses.R')
source('~/Documents/BIM/IntLim-master/IntLim/R/01_ReadData.R')
setwd("/Users/liumingrui/Documents/BIM/IntLim-master/IntLim/inst/extdata")


shinyUI(fluidPage(
    titlePanel("Uploading Files"),
    sidebarLayout(
        sidebarPanel(
            fileInput('file1', 'Choose CSV File',
                      accept=c('text/csv', 
                               'text/comma-separated-values,text/plain', 
                               '.csv')),
            textInput("metabid", "Metab ID", "BIOCHEMICAL"),
            textInput("geneid", "Gene ID", "X")
      
        ),
            
        mainPanel(
            pre(includeText("/Users/liumingrui/Documents/BIM/IntLim-master/IntLim/inst/extdata/textformat.txt")),
            pre(textOutput('contents'))
        )
    )
    )
)

