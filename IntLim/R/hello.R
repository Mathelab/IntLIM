# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

library(MultiAssayExperiment)
library(GenomicRanges)


path<-"/Users/liumingrui/Desktop/NCI60.sample.csv"
setwd("/Users/liumingrui/Documents/BIM/Files_for_understanding_Data_and_Package")


createIntlim<-function(path,log2=FALSE){
    
    stopifnot(is.character(path))
    
    if (!file.exists(path)) {
        stop("CSV input file does not exist")
    }
    
    csvfile <- as.data.frame(t(read.table(path, header=TRUE,
                                          sep=",", row.names="type")))
    
    if(toString(csvfile$geneData)==''){
        stop("GeneData is required")
    }
    MData<-as.data.frame(read.csv(file=toString(csvfile$metabData),row.names = 1))
    if(toString(csvfile$metabData)==''){
        stop("MetabData is required")
    }
    GData<-as.data.frame(read.csv(file=toString(csvfile$geneData),row.names = 1))
    if(toString(csvfile$pData)==''){
        stop("pData is required")
    }
    pData<-as.data.frame(read.csv(file=toString(csvfile$pData)))
    
    Mmetadata<-NULL
    Gmetadata<-NULL
    if(!(toString(csvfile$metabMetadata)=='')){
        Mmetadata<-as.data.frame(read.csv(file=toString(csvfile$metabMetadata)))
    }
    if(!toString(csvfile$geneMetadata) ==''){
        Gmetadata<-as.data.frame(read.csv(file=toString(csvfile$geneMetadata)))
    }
    
    ## ExperimentList
    MData<-t(meanmetab)
    GData <- genedat
    objlist <- list("Gene" = GData, "Metab" = MData)
    
    
    ## sampleMap
    Mmap <- data.frame(primary = colnames(MData), assay =colnames(MData))
    Gmap <- data.frame(primary = colnames(GData), assay = colnames(GData))
    listmap <- list(Gmap, Mmap)
    names(listmap) <- c("Gene", "Metab")
    dfmap <- listToMap(listmap)
    
    
    ##ColData
   
    colData<-unique(sampannot)
    rownames(colData)<-colData$sampannot.cell_line
    
    
    
    ##MetaData
    listMetaData<-list(Gmetadata,Mmetadata)
    names(listMetaData)<-c("Gene", "Metab")
    
    
    ## Create MultiAssayExperiment Object
    MultiAssay1 <- MultiAssayExperiment(objlist, colData, dfmap,listMetaData)
    return(MultiAssay1)
}



