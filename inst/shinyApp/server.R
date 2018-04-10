options(shiny.trace=F)
library(readr)
shinyServer(function(input, output, session) {

    #Desktop file input==================================================================================================
  fixUploadedFilesNames <- function(x) {
    if (is.null(x)) {
      return()
    }
    oldNames = x$datapath
    newNames = file.path(dirname(x$datapath),x$name)
    file.rename(from = oldNames, to = newNames)
    x$datapath <- newNames
    x
  }

  output$filename <- renderPrint({
    myFile <- fixUploadedFilesNames(input$file1) #fixing uploaded filenames.
    if (is.null(myFile)) {
      cat("Please name your input file as input.csv & select multiple files by clicking the browse button")
    } else if (length(myFile$name) !=6) {
      cat("please select upto 6 files")
    } else {
      paste(myFile$name,"- uploaded")
    }

  })

  output$idChooseM <- renderUI({

    if (is.null(input$file1$datapath)) {

    } else {
      myFile <- fixUploadedFilesNames(input$file1)
      indexfile = which(myFile$name == 'input.csv')
      if (indexfile == which(myFile$name == 'input.csv')) {
        file <- read_csv(myFile$datapath[indexfile])
        rownames(file) <- file$type
        if (file["metabMetaData", "filenames"] == "") {
          return()
        }
        textInput("metabid", "Metab ID", "")

      } else {

        tags$b(paste("Please rename input file as input.csv"))

      }
    }

  })

  output$idChooseG <- renderUI({

    if (is.null(input$file1)) {

    } else {
      myFile <- fixUploadedFilesNames(input$file1)
      indexfile = which(myFile$name == 'input.csv')
      if (indexfile == which(myFile$name == 'input.csv')) {
        file <- read_csv(myFile$datapath[indexfile])
        rownames(file) <- file$type
        if (file["metabMetaData", "filenames"] == "") {
          return()
        }
        textInput("geneid", "Gene ID", "")

      } else {

        tags$b(paste("Please rename input file as input.csv"))

      }
    }

  })

  ## Windows end
  multiData <- eventReactive(input$run, {
    print(Sys.time())
    myFile <- fixUploadedFilesNames(input$file1)
    print(length(myFile))
    indexfile = which(myFile$name == 'input.csv')
    IntLIM::ReadData(req(myFile$datapath[[indexfile]]),
                     input$metabid, input$geneid)

  })

  output$stats<-renderDataTable({
        table<- as.data.frame(t(IntLIM::ShowStats(multiData())))
        colnames(table)<-"value"
        cbind(names=rownames(table),table)

    },options = list(dom = 't'))


    output$plot<-renderUI(
        IntLIM::PlotDistributions(multiData())
    )

    #filter data==================================================================================================

    FmultiData<-eventReactive(input$run2,{
        if(input$run2==0){
            FmultiData<-multiData()
        }
        if(input$run2!=0){
            FmultiData<-IntLIM::FilterData(multiData(),
		geneperc=input$geneperc,
		metabperc=input$metabperc,
		metabmiss=input$metabmiss)
        }

        FmultiData
    },ignoreNULL=FALSE)

    output$Ostats<-renderDataTable({
        if(input$run2==0) return()
        table<- as.data.frame(t(IntLIM::ShowStats(multiData())))
        colnames(table)<-"value"
        cbind(names=rownames(table),table)

    },options = list(dom = 't'))


    output$Oplot<-renderUI({
        if(input$run2==0) return()
        IntLIM::PlotDistributions(multiData())
    }
    )
    output$Fstats<-renderDataTable({
        if(input$run2==0) return()
        table<- as.data.frame(t(IntLIM::ShowStats(FmultiData())))
        colnames(table)<-"value"
        cbind(names=rownames(table),table)

    },options = list(dom = 't'))


    output$Fplot<-renderUI({
        if(input$run2==0) return()
        IntLIM::PlotDistributions(FmultiData())
    }
    )

    output$downloadFdata <- downloadHandler(
        filename = "Filtered data.zip",
        content = function(con) {
            IntLIM::OutputData(FmultiData(),con)
        }
    )

    #run Lntlim==================================================================================================
    output$choosestype <- renderUI({

        choice<-reactive({
            Biobase::varLabels(FmultiData()[["expression"]])
        })

        selectInput("stype", "Sample Type:",
                    c(Choose='',choice()),selected = NULL)
    })


    myres <- eventReactive(input$run3,{
        shinyjs::html("text", "")
        IntLIM::RunIntLim(FmultiData(),stype=input$stype,outcome='metabolite')

    })
    diffcorr<-reactive(input$diffcorr1)
    pvalcutoff<-reactive(input$pvalcutoff1)

    output$volcanoPlot<-renderPlot(
        {IntLIM::pvalCorrVolcano(myres(),FmultiData(),input$nrpoints,diffcorr(),pvalcutoff())},
	height=500
    )
    output$Pdist<-renderPlot({

        IntLIM::DistPvalues(myres(),breaks = input$breaks)

    })

    output$Ptext<-renderPrint(

            if(!is.null(myres())){

        ("Distribution of unadjusted p-values (a peak close to zero suggests that there are significant gene:metabolite pairs that are found).")
            }

    )


    #heatmap==================================================================================================
    # observe({
    #     if(!is.null(input$diffcorr)&&!is.null(input$pvalcutoff)){
    #     diffcorr<-reactive(input$diffcorr)
    #     pvalcutoff<-reactive(input$pvalcutoff)
    #     }
    # })
    output$numericChoice1<-renderUI(
        numericInput("pvalcutoff","cutoff of FDR-adjusted p-value for filtering(0 - 1) :",pvalcutoff(), min = 0, max = 1)
    )
    output$numericChoice2<-renderUI(
        numericInput("diffcorr", "cutoff of differences in correlations for filtering (0-1):",diffcorr(), min = 0, max = 1)
    )

    myres2 <- eventReactive(input$run4,{

                IntLIM::ProcessResults(myres(),FmultiData(),pvalcutoff=pvalcutoff(),
                               diffcorr=diffcorr(),
                               corrtype=input$corrtype,
                               treecuts=input$treecuts)




    })

    output$heatmap<-plotly::renderPlotly({
        IntLIM::CorrHeatmap(myres2(),treecuts=input$treecuts)
    }
    )
    output$downloadData <- downloadHandler(
        filename = "results.csv",
        content = function(con) {
            IntLIM::OutputResults(req(myres2()),con)
        }
    )

    #scatter plot=============================================================================================

    pairTable<-reactive({
	mydat <- req(myres2())
	mydat@filt.results
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
    scatterrows<-eventReactive(input$run5,{
        input$table_rows_selected
    })
    #output$temp<-renderPrint(as.matrix(scatterrows()))

    output$scatterplot<-renderUI({
            a<-as.matrix(scatterrows())
            pair1<-as.matrix(pairTable()[a[1,],])
            geneName1<-pair1[,"gene"]
            metabName1<-pair1[,"metab"]
            splot1<-IntLIM::PlotGMPair(FmultiData(),input$stype,geneName=geneName1,metabName=metabName1)
            if(length(input$table_rows_selected) > 1){
                pair2<-as.matrix(pairTable()[a[2,],])
                geneName2<-pair2[,"gene"]
                metabName2<-pair2[,"metab"]
                splot2<-IntLIM::PlotGMPair(FmultiData(),input$stype,geneName=geneName2,metabName=metabName2)

                p <-htmltools::browsable(highcharter::hw_grid(splot1, splot2, ncol = 2, rowheight = 550))
            }
            else{
                    p<-htmltools::browsable(highcharter::hw_grid(splot1, ncol = 1, rowheight = 550))
            }
            return(p)
    })

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
                HTML(paste("Data is loaded.",
                           "You can proceed to Step 2 (optional) or Step 3.",
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
                HTML(paste("Data filtering is complete.",
                           "You can proceed to Step 3",
                           sep = "<br/>")),
                icon = icon("thumbs-up", lib = "glyphicon"),
                color = "green", fill = TRUE)
        }
    })

    output$statusbox3 <- renderInfoBox({
        if (!is.null(input$stype=="")) {
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
                HTML(paste("IntLIM models are calculated.",
                           "You can proceed to Step 4",
                           sep = "<br/>")),
                icon = icon("thumbs-up", lib = "glyphicon"),
                color = "green", fill = TRUE)
        }
    })

    output$statusbox4 <- renderInfoBox({
        if (input$run4==0) {
            infoBox(
                "Status",
                "Press Run button",
                icon = icon("flag", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}

        else if (!is.null(myres2())) {
            infoBox(
                "Status",
                HTML(paste("Results processed.",
                           "You can proceed to Step 5",
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
                "Press Run",
                icon = icon("flag", lib = "glyphicon"),
                color = "aqua",
                fill = TRUE
            )}

        else if (!is.null(scatterrows())) {
            infoBox(
                "Status",
                HTML(paste("Scatter plot running complete!",

                           sep = "<br/>")),
                icon = icon("thumbs-up", lib = "glyphicon"),
                color = "green", fill = TRUE)
        }


    })



})

