library(shiny)

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
        OutputStats(multiData())
    )
    output$plot<-renderPlot(
        PlotDistributions(multiData())
    )
    output$Fstats<-renderTable(
        OutputStats(FmultiData())
    )
    
    
    
})