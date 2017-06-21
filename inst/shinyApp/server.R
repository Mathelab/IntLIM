library(shiny)

shinyServer(function(input, output) {
    
    multiData <- reactive({
        inFile<-input$file1
        multiData<-ReadData(inFile$datapath,input$metabid,input$geneid)
        if(input$filter){
            multiData<-FilterData(multiData,geneperc=input$geneperc,metabperc=input$metabperc)
        }
        multiData
    })

    output$contents <- renderPrint({
        if (is.null(input$file1))
            return(NULL)
        multiData()
    })
    
    output$stats<-renderTable(
        
        OutputStats(multiData())
    )
    output$plot<-renderPlot(
        PlotDistributions(multiData())
      
    )


        
})
