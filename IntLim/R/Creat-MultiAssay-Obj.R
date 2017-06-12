library(MultiAssayExperiment)

#' Read in CSV file
#'
#' The metadata associated with data files to be analyzed in IntLim is supplied
#' as a CSV file. The software will automatically retrieve the file path of
#' input CSV so it is important that all analysis files are in the same folder
#' as CSV file. The input file need have standard format and the path for all the needed CSV file.
#'
#' @param csvPath csvPath
#' @param log2 boolean
#' @return MultiAssayObject of CSV file
#'
#' @examples
#' \dontrun{
#'  MyMultiAssay1<-creatIntlim("breast.sample.csv",log2=FALSE)
#'
#' @export
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
    pData<-as.data.frame(read.csv(file=toString(csvfile$pData),row.names = 1))

    Mmetadata<-NULL
    Gmetadata<-NULL
    if(!(toString(csvfile$metabMetadata)=='')){
        Mmetadata<-as.data.frame(read.csv(file=toString(csvfile$metabMetadata)))
    }
    if(!toString(csvfile$geneMetadata) ==''){
        Gmetadata<-as.data.frame(read.csv(file=toString(csvfile$geneMetadata)))
    }

    ## ExperimentList
    MData<-MData
    GData <- GData
    objlist <- list("Gene" = GData, "Metab" = MData)


    ## sampleMap
    Mmap <- data.frame(primary = colnames(MData), assay =colnames(MData))
    Gmap <- data.frame(primary = colnames(GData), assay = colnames(GData))
    listmap <- list(Gmap, Mmap)
    names(listmap) <- c("Gene", "Metab")
    dfmap <- listToMap(listmap)


    ##ColData

    colData<-unique(pData)
    rownames(colData)<-row.names(pData)



    ##MetaData
    listMetaData<-list(Gmetadata,Mmetadata)
    names(listMetaData)<-c("Gene", "Metab")


    ## Create MultiAssayExperiment Object
    MultiAssay1 <- MultiAssayExperiment(objlist, colData, dfmap,listMetaData)
    return(MultiAssay1)
}
}



