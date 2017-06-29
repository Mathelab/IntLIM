require(shiny)
require(Biobase)
require(highcharter)
require(magrittr)
require(rCharts)
require(IntLim)
options(shiny.trace=TRUE)

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
        if(input$dataset=="metabolite"){
        choice<-reactive({
            varLabels(FmultiData()[["metabolite"]])
         })
        }else if(input$dataset=="gene"){
            choice<-reactive({
            varLabels(FmultiData()[["expression"]])
            })
        }
        selectInput("stype", "Sample Type:", 
                    choices=choice())
    })
    

    myres <- reactive({
        RunIntLim(FmultiData(),stype=input$stype,outcome=input$dataset)
    })
    logText <- reactive({
        textoutput <- capture.output(data <- myres())
        
        
    })
    output$process<-renderPrint({
        logText()
        return(print(logText()))
    })
    output$Pdist<-renderPlot({
        DistPvalues(myres())
    })
    
    output$heatmap<-renderHighchart({
        myres2 <- IntLim::ProcessResults(myres(),FmultiData())
        CorrHeatmap(myres2)
    }
    )
    
    
    
})
