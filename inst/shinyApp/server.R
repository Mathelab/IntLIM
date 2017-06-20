library(shiny)

shinyServer(function(input, output) {
    
    multiData <- reactive({
        inFile<-input$file1
        multiData<-ReadData(inFile$datapath,input$metabid,input$geneid)
        multiData
        
    })
    output$contents <- renderPrint({
        if (is.null(input$file1))
            return(NULL)
        multiData()
        
    })
    output$group<-renderPrint({
        multiData<-multiData()
       
        if(input$groupA=="T"){
            groupA<-subset(multiData,, DIAG=="TUMOR")
        }else{
            groupA<-subset(multiData,, DIAG=="NORMAL")
        }
        
        if(input$groupB=="T"){
            groupB<-subset(multiData,, DIAG=="TUMOR")
        }else{
            groupB<-subset(multiData,, DIAG=="NORMAL")
        }
        
        A.gene.set<-groupA[["expression"]]
        A.metab.set<-groupA[["metabolite"]]
        B.gene.set<-groupB[["expression"]]
        B.metab.set<-groupB[["metabolite"]]
        
        A.gene.mat <- exprs(A.gene.set)
        B.gene.mat <- exprs(B.gene.set)
        A.metab.mat <- getMetabData(A.metab.set)
        B.metab.mat <- getMetabData(B.metab.set)
        
        p.val.bc.list <- get.p.vals.mlm(t(A.gene.mat), t(B.gene.mat), t(A.metab.mat), t(B.metab.mat))
        p.val.bc.list.adjust <- p.val.adjust(p.val.bc.list)
        rownames(p.val.bc.list.adjust) <- rownames(p.val.bc.list)
        colnames(p.val.bc.list.adjust) <- colnames(p.val.bc.list)
        
        p.val.bc.list.adjust
        
    })
        
        
})
