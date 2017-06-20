#' Filter data by abundance and/or CV
#'
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression, 
#' metabolite abundances, and associated meta-data
#' @param geneCV CV cutoff to be used for filtering genes (e.g. remove genes CV < 5)
#' @param metabCV CV cutoff to be used for filtering metabolites
#' @param geneperc percentile cutoff for filtering genes (e.g. remove genes with values < 15th
#'  percentile)
#' @param metabperc percentile cutoff for filtering metabolites
#' @return filtData MultiDataSet object with input data after filtering
#'
#' @examples
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "test.csv")
#' mydata <- ReadData(csvfile,metabid='BIOCHEMICAL',geneid='X')
#' mydatafilt <- FilterData(mydata,geneperc=15,metabCV=10)
#' @export
FilterData <- function(inputData,geneCV=NULL,metabCV=NULL,geneperc=NULL,
	metabperc=NULL) {

    # Check that input is a MultiDataSet
    if (class(inputData) != "MultiDataSet") {
	stop("input data is not a MultiDataSet class")
    }
    mytypes <- names(Biobase::assayData(inputData))
    if(!any(mytypes=="expression") || !any(mytypes=="metabolite")) {
	stop("input data must contain assayData of type 'metabolite' and 'expression.
	Try reading in the data with the ReadData function")
    }	

    # Check that at least one parameter is not null
    len <- length(c(geneCV,metabCV,geneperc,metabperc))
    if (len==0) {
        warning("All filtering parameters are NULL so the data remains unfiltered")
	return(inputData)
    }
    else {return(NULL)
	}

}

