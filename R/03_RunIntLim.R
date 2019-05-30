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
#' @param covar Additional variables from the phenotypic data that be integrated into linear model
#' @param class.covar Describing whether additional variables are 'numeric' or 'categorial'
#' @param continuous boolean to indicate whether the data is continuous or discrete
#' @return IntLimModel object with model results
#'
#' @examples
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' \dontrun{
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' }
#' @export
RunIntLim <- function(inputData,stype=NULL,outcome="metabolite", covar=NULL, class.covar=NULL, continuous = FALSE){


    if (class(inputData) != "MultiDataSet") {
        stop("input data is not a MultiDataSet class")
    }
    mytypes <- names(Biobase::assayData(inputData))
    if(!any(mytypes=="expression") || !any(mytypes=="metabolite")) {
        stop("input data must contain assayData of type 'metabolite' and 'expression.
        Try reading in the data with the ReadData function")
    }
    if (is.null(stype)) {stop("Please set the variable type (e.g. sample group)")}

#    incommon<-MultiDataSet::commonSamples(inputData)
#    mp <- as.character(Biobase::pData(incommon[["metabolite"]])[,stype])
#    gp <- as.character(Biobase::pData(incommon[["expression"]])[,stype])

    incommon <- getCommon(inputData,stype,covar,class.covar=class.covar)


    if(!continuous & length(unique(stats::na.omit(incommon$p))) != 2) {
	    stop(paste("IntLim currently requires only two categories.  Make sure the column",stype,"only has two unique values"))
    }

    print("Running the analysis on")
    print(table(incommon$p))

    ptm <- proc.time()

    myres <- RunLM(incommon,outcome=outcome,type=incommon$p,covar=covar, continuous = continuous)

    print(proc.time() - ptm)
    myres@stype=stype
    myres@outcome=outcome

    if(!is.null(covar)){
        covariate <- colnames(incommon$covar_matrix)
        class.var <- c()
        for(i in 1:length(covar)){
            class.var[i] <- class(incommon$covar_matrix[,i])
        }

        myres@covar <- data.frame(covariate,class.var)
    }
    return(myres)
}
