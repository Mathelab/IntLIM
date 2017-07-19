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
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' \dontrun{
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' }
#' @export
RunIntLim <- function(inputData,stype=NULL,outcome="metabolite"){

    if (class(inputData) != "MultiDataSet") {
        stop("input data is not a MultiDataSet class")
    }
    mytypes <- names(Biobase::assayData(inputData))
    if(!any(mytypes=="expression") || !any(mytypes=="metabolite")) {
        stop("input data must contain assayData of type 'metabolite' and 'expression.
        Try reading in the data with the ReadData function")
    }
    if (is.null(stype)) {stop("Please set the variable type (e.g. sample group)")}

    incommon<-MultiDataSet::commonSamples(inputData)
    mp <- as.character(Biobase::pData(incommon[["metabolite"]])[,stype])
    gp <- as.character(Biobase::pData(incommon[["expression"]])[,stype])
    if(all.equal(mp,gp)[1] != TRUE) {
	stop(paste("The column", stype,"for the samples in common between the metabolite and gene datasets are not equal.  Please check your input."))
    }
    mp[which(mp=="")]=NA
    if(length(unique(stats::na.omit(mp))) > 2) {
	stop(paste("IntLim currently requires only two categories.  Make sure the column",stype,"only has two unique values"))
    }
    # Get the samples where the outcome has a non-missing value
    if(any(is.na(mp))) {
	keepers <- which(!is.na(mp))
        namestokeep <- rownames(Biobase::pData(incommon[["metabolite"]]))[keepers]
    	incommon <- incommon[namestokeep]
	mp <- mp[keepers]
    }

    print("Running the analysis on")
    print(table(mp))

    ptm <- proc.time()
    myres <- RunLM(incommon,outcome=outcome,type=mp)
    print(proc.time() - ptm)
    myres@stype=stype
    myres@outcome=outcome
    return(myres)
}
