#' Filter data by abundance and/or CV
#'
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression, 
#' metabolite abundances, and associated meta-data
#' @param geneperc percentile cutoff for filtering genes (e.g. remove genes with mean values 
#' < 15th percentile)
#' @param metabperc percentile cutoff for filtering metabolites
#' @return filtData MultiDataSet object with input data after filtering
#'
#' @examples
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' inputData <- ReadData(csvfile,metabid='id',geneid='id')
#' inputDatafilt <- FilterData(inputData,geneperc=15)
#' @export
FilterData <- function(inputData,geneperc=NULL,metabperc=NULL) {

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
    len <- length(c(geneperc,metabperc))
    if (len==0) {
        warning("All filtering parameters are NULL so the data remains unfiltered")
	return(inputData)
    }
    else {
	mygenes <- Biobase::assayDataElement(inputData[["expression"]], 'exprs')
	mymetab <- Biobase::assayDataElement(inputData[["metabolite"]], 'metabData')
	if(!is.null(geneperc)) {
		if(geneperc>1) {geneperc=geneperc/100}
		mymean <- as.numeric(apply(mygenes,1, function(x)
			mean(x,na.rm=T)))
		keepers <- which(mymean > stats::quantile(mymean,geneperc))
		mygenes <- mygenes[keepers,]
		pgenes <- Biobase::pData(inputData[["expression"]])
		fgenes <- Biobase::fData(inputData[["expression"]])[keepers,]
	} else {
		mygenes <- mygenes
                fgenes <- Biobase::fData(inputData[["expression"]])
	}
	if(!is.null(metabperc)) {
		if(metabperc>1) {metabperc=metabperc/100}
                mymean <- as.numeric(apply(mymetab,1, function(x)
                        mean(x,na.rm=T)))
                keepers <- which(mymean > stats::quantile(mymean,metabperc))
                mymetab <- mymetab[keepers,]
		fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]])[keepers,])
        } else {
		mymetab <- mymetab
		fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]]))
	}
	# Now reconstruct the multidataset object
	gene.set <- Biobase::ExpressionSet(assayData=mygenes)
        Biobase::fData(gene.set) <- fgenes
        Biobase::pData(gene.set) <- Biobase::pData(inputData[["expression"]])
        
	metab.set <- methods::new("MetaboliteSet",metabData = mymetab,
                phenoData =  Biobase::AnnotatedDataFrame(data = Biobase::pData(inputData[["metabolite"]])), 
		featureData =  fmetab)

	multi <- MultiDataSet::createMultiDataSet()
        multi1 <- MultiDataSet::add_genexp(multi, gene.set)
        filtdata <- add_metabolite(multi1, metab.set)

	return(filtdata)
    }
}

