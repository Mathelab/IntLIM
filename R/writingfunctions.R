#' Output data into individual CSV files.  All data will be zipped into one file with all data.
#'
#' @param inputData data output from ReadData() or FilterData() function
#' @param filename name of file to be output (default: '~/output.zip') 
#'
#' @return the filename of the CSV file with results named with cohort
#'
#' @examples
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' inputData <- ReadData(csvfile,metabid='id',geneid='id')
#' inputDatafilt <- FilterData(inputData,geneperc=0.5)
#' OutputData(inputData=inputDatafilt,filename="~/FilteredData.zip")
#' @export

OutputData <- function (inputData=NULL,filename="~/output.zip"){
        mygenes <- Biobase::assayDataElement(inputData[["expression"]], 'exprs')
        mymetab <- Biobase::assayDataElement(inputData[["metabolite"]], 'metabData')
        fgenes <- Biobase::fData(inputData[["expression"]])
	fmetab <- Biobase::fData(inputData[["metabolite"]])
	phenoData <- Biobase::pData(inputData[["metabolite"]])
	utils::write.csv(mygenes,"~/GeneData.csv",quote=F)
        utils::write.csv(mymetab,"~/MetabData.csv",quote=F)
        utils::write.csv(fgenes,"~/MetaGenes.csv",quote=F)
        utils::write.csv(fmetab,"~/MetaMetab.csv",quote=F)
        utils::write.csv(phenoData,"~/MetaSamples.csv",quote=F)
	utils::zip(zipfile=path.expand(filename),
		files=c(path.expand("~/GeneData.csv"),
			path.expand("~/MetabData.csv"),
			path.expand("~/MetaGenes.csv"),
			path.expand("~/MetaMetab.csv"),
			path.expand("~/MetaSamples.csv")))
	file.remove(path.expand("~/GeneData.csv"))
	file.remove(path.expand("~/MetabData.csv"))
	file.remove(path.expand("~/MetaGenes.csv"))
	file.remove(path.expand("~/MetaMetab.csv"))
	file.remove(path.expand("~/MetaSamples.csv"))
}



#' Output results into a zipped CSV file.  Results include gene and metabolite pairs, along with model interaction p-values, and correlations in each group being evaluated.
#'
#' @param inputResults IntLimResults object with model results (output of ProcessResults())
#' @param filename name of file to be output (default: '~/results.csv')
#'
#' @return the filename of the CSV file with results named with cohort
#'
#' @examples
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' inputData <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(inputData,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,inputData)
#' OutputResults(inputResults=myres,filename="~/results.csv")
#' @export

OutputResults <- function (inputResults=NULL,filename="~/results.csv"){
	if(is.null(inputResults)) {stop("Input results from ProcessResults()")}
	if(is.null(inputResults@filt.results)) {stop("You must run ProcessResults() on results table first")}
	utils::write.csv(inputResults@filt.results,filename,quote=T,row.names=F)
}


