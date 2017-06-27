library(shiny)
require(shinydashboard)
require(magrittr)

require(rCharts)
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
    })
    
   
    output$stats<-renderTable(
        as.matrix(OutputStats(multiData()))
    )
    
    output$plot<-renderHighchart(
        PlotDistributions(multiData())
       )
    

    
    
    output$Fstats<-renderTable(
        OutputStats(FmultiData())
    )
    output$Fplot<-renderPrint(
        PlotDistributions(FmultiData())
        
    )
    
    output$Pdist<-renderPlot({
        myres <- RunIntLim(FmultiData(),stype="DIAG")
        DistPvalues(myres)
    })
    
    output$heatmap<-renderPrint({
        myres2 <- ProcessResults(myres,inputData)
        CorrHeatmap(myres)
        }
    )
    
    
    
})