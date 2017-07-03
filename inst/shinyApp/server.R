
options(shiny.trace=TRUE)

shinyServer(function(input, output) {
    
    # based on the Shiny fileInput function
    fileInput2 <- function(inputId, label = NULL, labelIcon = NULL, multiple = FALSE,accept = NULL, width = NULL, progress = TRUE, ...) {
        # add class fileinput_2 defined in UI to hide the inputTag
        inputTag <- tags$input(id = inputId, name = inputId, type = "file", 
                               class = "fileinput_2")
        if (multiple) 
            inputTag$attribs$multiple <- "multiple"
        if (length(accept) > 0) 
            inputTag$attribs$accept <- paste(accept, collapse = ",")
        
        div(..., style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), ";"), 
            inputTag,
            # label customized with an action button
            tags$label(`for` = inputId, div(icon(labelIcon), label, 
                                            class = "btn btn-default action-button")),
            # optionally display a progress bar
            if(progress)
                tags$div(id = paste(inputId, "_progress", sep = ""), 
                         class = "progress shiny-file-input-progress", 
                         tags$div(class = "progress-bar")
                )
        )
    }          
    #file input
    output$filename<-renderPrint(
        {
            inFile <- input$file1
            if (is.null(inFile)){
                cat("Please chose file")
            }else{
            paste("File name:", inFile$name)
            }
        }
    )
    multiData <- eventReactive(input$run,{
        
        inFile<-input$file1
        multiData<-ReadData(inFile$datapath,input$metabid,input$geneid)
        multiData
    })
    
    output$stats<-renderTable(
        t(OutputStats(multiData())),
        include.rownames=TRUE,
        include.colnames=FALSE
    )
    
    output$plot<-renderUI(
        PlotDistributions(multiData())
    )
    #filter data
    FmultiData<-eventReactive(input$run2,{
        FmultiData<-multiData()
        if(input$filter){
            FmultiData<-FilterData(multiData(),geneperc=input$geneperc,metabperc=input$metabperc)
        }
        FmultiData
    },ignoreNULL=FALSE)
    
    output$Fstats<-renderTable(
        as.data.frame(t(OutputStats(FmultiData()))),
        include.rownames=TRUE,
        include.colnames=FALSE
    )
    output$Fplot<-renderUI(
        PlotDistributions(FmultiData())
        
    )
    
    
    #adjusted p values
    output$choosestype <- renderUI({
        if(input$dataset=="metabolite"){
            choice<-reactive({
                varLabels(FmultiData()[[input$dataset]])
            })
        }else if(input$dataset=="gene"){
            choice<-reactive({
                varLabels(FmultiData()[["expression"]])
            })
        }
        selectInput("stype", "Sample Type:", 
                    choices=choice())
    })
        
       
    myres <- eventReactive(input$run3,{
       RunIntLim(FmultiData(),stype=input$stype,outcome=input$dataset)
    })
    output$Pdist<-renderPlot({
       DistPvalues(myres()@interaction.adj.pvalues)
        
    })
   
    
    #heatmap
    
    myres2 <- eventReactive(input$run4,{
        IntLim::ProcessResults(myres(),FmultiData())
    })
    output$heatmap<-renderHighchart({
        CorrHeatmap(myres2())
    }
    )
    
    
    
})
