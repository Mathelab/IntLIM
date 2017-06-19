library(shiny)

shinyServer(function(input, output) {
    

    output$contents <- renderPrint({
        
        
        
        inFile <- input$file1
        if (is.null(inFile))
            return(NULL)
        
        multiData<-ReadData(inFile$datapath,input$metabid,input$geneid)
        multiData
    })
    
        
        
})
