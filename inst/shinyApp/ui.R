
library(shiny)
library(MultiAssayExperiment)
setwd("/Users/liumingrui/Documents/BIM/IntLiMpackage")
source('~/Documents/BIM/IntLim-master/R/MetaboliteSet_addMetabolite.R')
source('~/Documents/BIM/IntLim-master/R/internalfunctions.R')
source('~/Documents/BIM/IntLim-master/R/01_ReadData.R')
source('~/Documents/BIM/IntLim-master/R/00_NewClasses.R')

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

