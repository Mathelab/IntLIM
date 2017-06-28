library(shiny)

require(Biobase)
require(highcharter)
require(magrittr)
require(rCharts)
require(IntLim)
shinyServer(function(input, output) {
    
    multiData <- eventReactive(input$run,{
        
        inFile<-input$file1
        multiData<-ReadData(inFile$datapath,input$metabid,input$geneid)
        multiData
    })
    
    
    FmultiData<-eventReactive(input$run2,{
        FmultiData<-multiData()
        if(input$filter){
            FmultiData<-FilterData(multiData(),geneperc=input$geneperc,metabperc=input$metabperc)
        }
        FmultiData
    },ignoreNULL=FALSE)
    
    
    output$stats<-renderTable(
        OutputStats(multiData())
    )
    
    output$plot<-renderUI(
        PlotDistributions(multiData())
    )
    
    
    
    
    output$Fstats<-renderTable(
        OutputStats(FmultiData())
    )
    output$Fplot<-renderUI(
        PlotDistributions(FmultiData())
        
    )
    
    
    
    output$choosestype <- renderUI({
        choice<-reactive({
            varLabels(FmultiData()[[input$dataset]])
        })
        selectInput("stype", "Sample Type:", 
                    choices=choice())
    })
    

    myres <- reactive({
        RunIntLim(FmultiData(),stype=input$stype,outcome=input$dataset)
    })
    
    output$Pdist<-renderPlot({
        DistPvalues(myres())
    })
    
    output$heatmap<-renderHighchart({
        myres2 <- IntLim::ProcessResults(myres(),FmultiData())
        plotCorrHeatmap(myres2)
    }
    )
    
    
    
})
