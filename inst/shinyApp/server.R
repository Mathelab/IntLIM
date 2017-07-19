
options(shiny.trace=TRUE)

shinyServer(function(input, output) {      
    #file input
    output$filename<-renderPrint(
        {
            inFile <- input$file1
            if (is.null(inFile)){
                cat("Please select CSV file by clicking the button above")
            }else{
            paste("File name:", inFile$name)
            }
        }
    )
    
    
    multiData <- eventReactive(input$run,{
        
        inFile<-input$file1
     
        multiData<-IntLim::ReadData(inFile$datapath,input$metabid,input$geneid)
        multiData
    })
    
    
    output$stats<-renderTable(
        t(IntLim::OutputStats(multiData())),
        include.rownames=TRUE,
        include.colnames=FALSE
    )
    
    output$plot<-renderUI(
        IntLim::PlotDistributions(multiData())
    )
    
    #filter data
    FmultiData<-eventReactive(input$run2,{
        FmultiData<-multiData()
        if(input$filter){
            FmultiData<-IntLim::FilterData(multiData(),geneperc=input$geneperc,metabperc=input$metabperc)
        }
        FmultiData
    },ignoreNULL=FALSE)
    
    output$Fstats<-renderTable(
        as.data.frame(t(IntLim::OutputStats(FmultiData()))),
        include.rownames=TRUE,
        include.colnames=FALSE
    )
    output$Fplot<-renderUI(
        IntLim::PlotDistributions(FmultiData())
        
    )
    
    
    #adjusted p values
    output$choosestype <- renderUI({
        
            choice<-reactive({
                Biobase::varLabels(FmultiData()[["expression"]])
            })
        
        selectInput("stype", "Sample Type:", 
                    choices=choice())
    })
        
       
    myres <- eventReactive(input$run3,{
        shinyjs::html("text", "")
        IntLim::RunIntLim(FmultiData(),stype=input$stype,outcome=input$dataset)
        
    })
    output$Pdist<-renderPlot({
        
        IntLim::DistPvalues(myres()@interaction.adj.pvalues)
        
    })
   
    
    #heatmap
    
    myres2 <- eventReactive(input$run4,{
        IntLim::ProcessResults(myres(),FmultiData())
    })
    output$heatmap<-highcharter::renderHighchart({
        IntLim::CorrHeatmap(myres2())
    }
    )
    
    #scatter plot
    output$table<-renderTable({
        a<-myres2()@corr
        b<-abs(a[,3]-a[,4])
        a$diff<-b
        a<-a[,-3]
        a<-a[,-3]
        top.type<-head(a[order(a[,3],decreasing = TRUE),])
        return(top.type)
    })
    
    
    
    output$chooseMetabName <- renderUI({
        
        metabName<-reactive({
            a<-myres2()@corr
            a[,1]
        })
        
        selectInput('metabName', 'Metab Name', choices=metabName(), selectize=TRUE)
    })
    
    output$chooseGeneName <- renderUI({
        
        geneName<-reactive({
            a<-myres2()@corr
            a[,2]
        })
        
        selectInput('geneName', 'Gene Name', choices=geneName(), selectize=TRUE)
    })
    
    stypeList<-eventReactive(input$run5,{
        s2<-input$stype
        expression<-FmultiData()[["expression"]]
        expression[[s2]]
        
    })
    
  
   
    output$scatterPlot<-highcharter::renderHighchart({
        IntLim::PlotGMPair(FmultiData(),stypeList(),geneName=input$geneName,metabName=input$metabName)
        })
    
    
    
})
