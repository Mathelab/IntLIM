
library(shiny)
library(MultiDataSet)
source('~/Documents/BIM/IntLim-master/IntLim/R/MetaboliteSet_addMetabolite.R')
source('~/Documents/BIM/IntLim-master/IntLim/R/internalfunctions.R')
source('~/Documents/BIM/IntLim-master/IntLim/R/AllClasses.R')
source('~/Documents/BIM/IntLim-master/IntLim/R/01_ReadData.R')

source('~/Documents/BIM/IntLiMpackage/IntLiMModelcode.R')
source('/Users/liumingrui/Documents/BIM/Files_for_understanding_Data_and_Package/GEMfunctionsv3.R')

setwd("/Users/liumingrui/Documents/BIM/IntLim-master/IntLim/inst/extdata")


shinyUI(fluidPage(
    
    titlePanel("Analyze metablon"),
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
            
        ),
            
        mainPanel(
            tabsetPanel(
                
                tabPanel('Data summary',pre(includeText("/Users/liumingrui/Documents/BIM/IntLim-master/IntLim/inst/extdata/textformat.txt"),
                                            pre(textOutput('contents'))
                                            )),
                tabPanel('Analyze By Group',pre(textOutput('group')))
                
            )
        )
    )
    )
)

