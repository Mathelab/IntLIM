options(shiny.trace=F)


shinyServer(function(input, output, session) {      
    #file input==================================================================================================

   rootVolumes <- c(Home = normalizePath("~"), getVolumes()(), WD = '.')
 #    rootVolumes <- c(Home = "/Users/ewymathe/Downloads/IntLim/vignettes/NCI-60", getVolumes()(), WD = '.')   

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

	myFile <- as.character(
             parseFilePaths(
               rootVolumes,
               input$file1)$datapath)
	if (is.null(myFile)){
                cat("Please select CSV file by clicking the button above")
            }else{
            cat(paste0("File to load:", myFile))
            }

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
    
    
    
    multiData <- eventReactive(input$run,{
        
        IntLim::ReadData(req(as.character(
            parseFilePaths(
                rootVolumes,
                input$file1)$datapath)),
            input$metabid,input$geneid)
       
        
    })
    
   
    
    output$stats<-renderDataTable({
        
        table<- as.data.frame(t(IntLim::OutputStats(multiData())))
        colnames(table)<-"value"
        cbind(names=rownames(table),table)
        
    },options = list(dom = 't'))
    
    
    output$plot<-renderUI(
        IntLim::PlotDistributions(multiData())
    )
    
    #filter data==================================================================================================
  
    FmultiData<-eventReactive(input$run2,{
        if(input$run2==0){
            FmultiData<-multiData()
        }
        if(input$run2!=0){
            FmultiData<-IntLim::FilterData(multiData(),geneperc=input$geneperc,metabperc=input$metabperc,metabmiss=input$metabmiss)
        }
        
        FmultiData
    },ignoreNULL=FALSE)
    
    output$Ostats<-renderDataTable({
        if(input$run2==0) return()
        table<- as.data.frame(t(IntLim::OutputStats(multiData())))
        colnames(table)<-"value"
        cbind(names=rownames(table),table)
        
    },options = list(dom = 't'))
    
    
    output$Oplot<-renderUI({
        if(input$run2==0) return()
        IntLim::PlotDistributions(multiData())
    }
    )
    output$Fstats<-renderDataTable({
        if(input$run2==0) return()
        table<- as.data.frame(t(IntLim::OutputStats(FmultiData())))
        colnames(table)<-"value"
        cbind(names=rownames(table),table)
        
    },options = list(dom = 't'))
    
    
    output$Fplot<-renderUI({
        if(input$run2==0) return()
        IntLim::PlotDistributions(FmultiData())
    }
    )
    
    
    #adjusted p values==================================================================================================
    output$choosestype <- renderUI({
        
        choice<-reactive({
            Biobase::varLabels(FmultiData()[["expression"]])
        })
        
        selectInput("stype", "Sample Type:", 
                    c(Choose='',choice()),selected = NULL)
    })
    
    
    myres <- eventReactive(input$run3,{
        shinyjs::html("text", "")
        IntLim::RunIntLim(FmultiData(),stype=input$stype,outcome=input$dataset)
        
    })
    output$Pdist<-renderPlot({
        
        IntLim::DistPvalues(myres())
        
    })
    
    
    #heatmap==================================================================================================
    
    myres2 <- eventReactive(input$run4,{
        IntLim::ProcessResults(myres(),FmultiData())
    })
    output$heatmap<-plotly::renderPlotly({
        IntLim::CorrHeatmap(myres2())
    }
    )
    
    #scatter plot=============================================================================================
    
    pairTable<-reactive({
        a<-myres2()@corr
        b<-round(abs(a[,3]-a[,4]),3)
        a$diff<-b
        a$PBO<-round(a$PBO,2)
        a$Leukemia<-round(a$Leukemia,2)
       
        table<-a[order(a[,5],decreasing = TRUE),]
        table
    })
    # reset <- reactiveValues(sel = "")
    # output$table<-DT::renderDataTable({
    #     input$table_rows_selected
    #     DT::datatable(as.matrix(pairTable()),selection = list(mode = 'multiple', selected = reset$sel))
    #     observe({
    #         if(length(input$table_rows_selected) > 2){
    #             reset$sel <- setdiff(input$table_rows_selected, input$table_row_last_clicked)
    #         }else{
    #             reset$sel <- input$table_rows_selected
    #         }
    #     })
    # 
    #     })
    
    output$table<-DT::renderDataTable(
        pairTable()
    )
    rows<-eventReactive(input$run5,{
        input$table_rows_selected
    })
    output$temp<-renderPrint(as.matrix(rows()))
    
    output$scatterplot<-renderUI({
        
            a<-as.matrix(rows())
            pair1<-as.matrix(pairTable()[a[1,],])
            geneName1<-pair1[,"gene"]
            metabName1<-pair1[,"metab"]
            splot1<-IntLim::PlotGMPair(FmultiData(),input$stype,geneName=geneName1,metabName=metabName1)
            
            
                if(length(input$table_rows_selected) > 1){
                    pair2<-as.matrix(pairTable()[a[2,],])
                    geneName2<-pair2[,"gene"]
                    metabName2<-pair2[,"metab"]
                    splot2<-IntLim::PlotGMPair(FmultiData(),input$stype,geneName=geneName2,metabName=metabName2)
            
                     p <-htmltools::browsable(highcharter::hw_grid(splot1, splot2, ncol = 2, rowheight = 550))
                }
                else{
                    p<-htmltools::browsable(highcharter::hw_grid(splot1, ncol = 1, rowheight = 550))
                }
                
        
            return(p)    
        
            
            
        
    })
    # output$scatterPlot1<-highcharter::renderHighchart({
    #     a<-as.matrix(rows())
    #     pair<-as.matrix(pairTable()[a[1,],])
    #     geneName<-pair[,"gene"]
    #     metabName<-pair[,"metab"]
    #     IntLim::PlotGMPair(FmultiData(),input$stype,geneName=geneName,metabName=metabName)
    #     
    # })
    # output$scatterPlot2<-highcharter::renderHighchart({
    #     a<-as.matrix(rows())
    #     if(length(rows())>1){
    #     pair<-as.matrix(pairTable()[a[2,],])
    #     geneName<-pair[,"gene"]
    #     metabName<-pair[,"metab"]
    #     IntLim::PlotGMPair(FmultiData(),input$stype,geneName=geneName,metabName=metabName)
    #     }else{
    #         highcharter::chart.renderer.text('You could see the second scatter plot by clicking another row')
    #     }
    #     
    #     
    # })
    
    #infobox
    output$statusbox1 <- renderInfoBox({
        if (is.null(input$file1)) {
            infoBox(
                "Status",
                "File Not Loaded Yet!",
                icon = icon("import", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}
        else if (!is.null(input$file1)&&input$run==0) {
            infoBox(
                "Status",
                "Step 1 is Not Complete Yet!",

                "Press Run button",

                icon = icon("warning-sign", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}
        else if (!input$run==0) {
            infoBox(
                "Status",
                HTML(paste("Object analyze Complete.",
                           "You can proceed to step2 (optional) or proceed to step3.",
                           sep = "<br/>")),
                icon = icon("thumbs-up", lib = "glyphicon"),
                color = "green", fill = TRUE)
        }
    })
    
    output$statusbox2 <- renderInfoBox({
        if (input$geneperc==0&&input$metabperc==0) {
            infoBox(
                "Status",
                "Please provide input",
                icon = icon("flag", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}
        else if (input$run2==0) {
            infoBox(
                "Status",
                "Press Run button",
                icon = icon("flag", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}
        else if (!input$run2==0) {
            infoBox(
                "Status",
                HTML(paste("Data filter complete.",
                           "You can proceed to step3",
                           sep = "<br/>")),
                icon = icon("thumbs-up", lib = "glyphicon"),
                color = "green", fill = TRUE)
        }
    })
    
    output$statusbox3 <- renderInfoBox({
        if (input$stype=="") {
            infoBox(
                "Status",
                "Please select your sample type",
                
                icon = icon("flag", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}
        else if (input$run3==0) {
            infoBox(
                "Status",
                "Step 3 is Not Complete Yet!",

               

                "Press Run button",
                "This function can take several minutes, please be patient",
                icon = icon("warning-sign", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}
        else if (!is.null(myres())) {
            infoBox(
                "Status",
                HTML(paste("IntLim models are calculated.",
                           "You can proceed to step4",
                           sep = "<br/>")),
                icon = icon("thumbs-up", lib = "glyphicon"),
                color = "green", fill = TRUE)
        }
    })
    
    output$statusbox4 <- renderInfoBox({
        if (input$run4==0) {
            infoBox(
                "Status",
                "Press Run to see the heatmap",
                icon = icon("flag", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}
        
        else if (!is.null(myres2())) {
            infoBox(
                "Status",
                HTML(paste("Heatmap running complete.",
                           "You can proceed to step5",
                           sep = "<br/>")),
                icon = icon("thumbs-up", lib = "glyphicon"),
                color = "green", fill = TRUE)
        }
        
        
    })
    output$statusbox5 <- renderInfoBox({
        if (is.null(input$table_rows_selected)) {
            infoBox(
                "Status",
                "Please choose a pair of gene and metab by clicking the table",
                icon = icon("flag", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}
        else if (input$run5==0) {
            infoBox(
                "Status",
                "Please Press Run button to run scatter plot",
                icon = icon("flag", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}
        
        else if (!is.null(rows())) {
            infoBox(
                "Status",
                HTML(paste("Scatter plot running complete!",
                           
                           sep = "<br/>")),
                icon = icon("thumbs-up", lib = "glyphicon"),
                color = "green", fill = TRUE)
        }
        
        
    })
    
    
    
})

