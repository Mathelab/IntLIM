options(shiny.trace=FALSE)


shinyServer(function(input, output, session) {      
    
    rootVolumes <- c(Home = normalizePath("~"), getVolumes()(), WD = '.')
    
    shinyFileChoose(input,'file1',
                    roots = rootVolumes,
                    session = session)
    
    output$filename <- renderPrint( {
        myFile <- as.character(
            parseFilePaths(
                rootVolumes,
                input$file1)$datapath)
        if (is.null(myFile)){
            cat("Please select CSV file by clicking the button above")
        }else{
            paste("File name:", myFile)
        }
    }
    )
    
    
    output$idChooseM <- renderUI({
        if (is.null(input$file1)){
            
        }else{
            myFile <- as.character(
                parseFilePaths(
                    rootVolumes,
                    input$file1)$datapath)
            file<- read.csv(myFile)
            rownames(file)<-file$type
            if(file["metabMetaData","filenames"]==""){
                return()
            }
            textInput("metabid", "Metab ID", "")
        }
        
        
    })
    
    output$idChooseG <- renderUI({
        if (is.null(input$file1)){
            
        }else{
            myFile <- as.character(
                parseFilePaths(
                    rootVolumes,
                    input$file1)$datapath)
            file<- read.csv(myFile)
            rownames(file)<-file$type
            if(file["geneMetaData","filenames"]==""){
                return()
            }
            textInput("geneid", "Gene ID", "")
        }
        
        
    })
    
    
    
    #file input
    #    output$filename<-renderPrint(
    #        {
    #            inFile <- xinput$file1
    #            if (is.null(inFile)){
    #                cat("Please select CSV file by clicking the button above")
    #            }else{
    #            paste("File name:", inFile$name)
    #            }
    #        }
    #    )
    
    
    multiData <- eventReactive(input$run,{
        
        multiData<-IntLim::ReadData(req(as.character(
            parseFilePaths(
                rootVolumes,
                input$file1)$datapath)),
            input$metabid,input$geneid)
        IntLim::FilterData(multiData,geneperc=5,metabperc=5)
        
    })
    
    
    
    output$stats<-renderDataTable({
        
        table<- as.data.frame(t(IntLim::OutputStats(multiData())))
        colnames(table)<-"value"
        cbind(names=rownames(table),table)
        
    },options = list(dom = 'ft'))
    
    
    output$plot<-renderUI(
        IntLim::PlotDistributions(multiData())
    )
    
    #filter data
    
    temp<-reactive(
        IntLim::FilterData(multiData(),geneperc=input$geneperc,metabperc=input$metabperc)
    )
    output$temp2<-renderPrint(temp())
    FmultiData<-eventReactive(input$run2,{
        if(input$run2==0){
            FmultiData<-multiData()
        }
        if(input$run2!=0){
            FmultiData<-IntLim::FilterData(multiData(),geneperc=input$geneperc,metabperc=input$metabperc)
        }
        
        FmultiData
    },ignoreNULL=FALSE)
    
    output$Fstats<-renderDataTable({
        table<- as.data.frame(t(IntLim::OutputStats(FmultiData())))
        colnames(table)<-"value"
        cbind(names=rownames(table),table)
        
    },options = list(dom = 'ft'))
    
    
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