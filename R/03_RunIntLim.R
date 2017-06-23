#' Run linear models and retrieve relevant statistics
#'
#' @include internalfunctions.R
#'
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param stype column name that represents sample type (by default, it will be used 
#' in the interaction term). Only 2 categories are currently supported.
#' @param outcome 'metabolite' or 'gene' must be set as outcome/independent variable
#' (default is 'metabolite')
#' @return IntLimModel object with model results
#'
#' @examples
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "test.csv")
#' mydata <- ReadData(csvfile,metabid='BIOCHEMICAL',geneid='X')
#' myres <- RunIntLim(mydata,stype="DIAG")
#' @export
ReadData <- function(inputFile,stype=NULL,outcome="metabolite"){

    if (class(inputData) != "MultiDataSet") {
        stop("input data is not a MultiDataSet class")
    }
    mytypes <- names(Biobase::assayData(inputData))
    if(!any(mytypes=="expression") || !any(mytypes=="metabolite")) {
        stop("input data must contain assayData of type 'metabolite' and 'expression.
        Try reading in the data with the ReadData function")
    }
    if (is.null(type)) {stop("Please set the variable type (e.g. sample group)")}

    incommon<-MultiDataSet::commonSamples(inputData)
    mp <- pData(incommon[["metabolite"]])[,stype]
    gp <- pData(incommon[["expression"]])[,stype]
    if(all.equal(mp,gp)[1] != TRUE) {
	stop(paste("The column", stype,"for the samples in common between the metabolite and gene datasets are not equal.  Please check your input."))
    }

    myres <- RunLM(incommon,outcome=outcome,type=mp)
}
