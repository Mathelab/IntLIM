library(shiny)

require(highcharter)
require(magrittr)
require(rCharts)
require(IntLim)
shinyServer(function(input, output) {
    
    multiData <- eventReactive(input$run,{
        
        inFile<-input$file1
        multiData<-IntLim::ReadData(inFile$datapath,input$metabid,input$geneid)
        multiData
    })
    
    
    FmultiData<-eventReactive(input$run2,{
        FmultiData<-multiData()
        if(input$filter){
            FmultiData<-IntLim::FilterData(multiData(),geneperc=input$geneperc,metabperc=input$metabperc)
        }
        FmultiData
    },ignoreNULL=FALSE)
    
    
    output$stats<-renderTable(
        as.matrix(IntLim::OutputStats(multiData()))
    )
    
    output$plot<-renderHighchart2(
        IntLim::PlotDistributions(multiData())
    )
    
    
    
    
    output$Fstats<-renderTable(
        IntLim::OutputStats(FmultiData())
    )
    output$Fplot<-renderHighchart(
        IntLim::PlotDistributions(FmultiData())
        
    )
    
    
    
    output$choosestype <- renderUI({
        choice<-reactive({
            Biobase::varLabels(FmultiData()[[input$dataset]])
        })
        selectInput("stype", "Sample Type:", 
                    choices=choice())
    })
    

    myres <- reactive({
        IntLim::RunIntLim(FmultiData(),stype=input$stype,outcome=input$dataset)
    })
    
    output$Pdist<-renderPlot({
        IntLim::DistPvalues(myres())
    })
    
    output$heatmap<-renderHighchart({
        myres2 <- IntLim::ProcessResults(myres(),FmultiData())
        IntLim::CorrHeatmap(myres2)
    }
    )
    
    
    
})
