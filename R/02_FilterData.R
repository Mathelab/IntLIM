#' Filter input data by abundance values (gene and metabolite data) and number of missing values (metabolite data only).
#'
#' Filter data by abundance (with user-input percentile cutoff) of missing values (with user-input percent cutoff). Missing values are commonly found in metabolomics data so the parameter currently only applies to metabolomics data.
#'
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression, 
#' metabolite abundances, and associated meta-data
#' @param geneperc percentile cutoff (0-1) for filtering genes (e.g. remove genes with mean values 
#' < 'geneperc' percentile) (default: 0)
#' @param metabperc percentile cutoff (0-1) for filtering metabolites (default: no filtering of metabolites) (default:0)
#' @param metabmiss missing value percent cutoff (0-1) for filtering metabolites (metabolites with > 80\% missing values will be removed) (default:0)
#' @return filtData MultiDataSet object with input data after filtering
#'
#' @examples
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' inputData <- ReadData(csvfile,metabid='id',geneid='id')
#' inputDatafilt <- FilterData(inputData,geneperc=0.5)
#' @export
FilterData <- function(inputData,geneperc=0,metabperc=0, metabmiss=0) {

    # Check that input is a MultiDataSet
    if (class(inputData) != "MultiDataSet") {
	stop("input data is not a MultiDataSet class")
    }
    mytypes <- names(Biobase::assayData(inputData))
    if(!any(mytypes=="expression") || !any(mytypes=="metabolite")) {
	stop("input data must contain assayData of type 'metabolite' and 'expression.
	Try reading in the data with the ReadData function")
    }	

    if(!is.null(geneperc) && geneperc > 1) {stop("geneperc parameter must be between 0 and 1")}
    if(!is.null(metabperc) && metabperc > 1) {stop("metabperc parameter must be between 0 and 1")}
    if(!is.null(metabmiss) && metabmiss > 1) {stop("metabmiss parameter must be between 0 and 1")}


    # Check that at least one parameter is not null
    len <- length(c(geneperc,metabperc,metabmiss))
    if ((geneperc+metabperc+metabmiss) ==0) {
        warning("All filtering parameters are NULL so the data remains unfiltered")
	return(inputData)
    }
    else {
	mygenes <- Biobase::assayDataElement(inputData[["expression"]], 'exprs')
	mymetab <- Biobase::assayDataElement(inputData[["metabolite"]], 'metabData')
	if(geneperc > 0) {
		if(geneperc>1) {geneperc=geneperc}
		mymean <- as.numeric(apply(mygenes,1, function(x)
			mean(x,na.rm=T)))
		keepers <- which(mymean > stats::quantile(mymean,geneperc))
		mygenes <- mygenes[keepers,]
		pgenes <- Biobase::pData(inputData[["expression"]])
		fgenes <- Biobase::fData(inputData[["expression"]])[keepers,]
	} else {
                print("No gene filtering is applied")
		mygenes <- mygenes
                fgenes <- Biobase::fData(inputData[["expression"]])
	}
	if(metabperc > 0) {
		if(metabperc>1) {metabperc=metabperc}
                mymean <- as.numeric(apply(mymetab,1, function(x)
                        mean(x,na.rm=T)))
                keepers <- which(mymean > stats::quantile(mymean,metabperc))
                mymetab <- mymetab[keepers,]
		fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]]))[keepers,]
        } else {
                print("No metabolite filtering by percentile is applied")
		mymetab <- mymetab
		fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]]))
	}
	if(metabmiss > 0) {
		missnum <- as.numeric(apply(mymetab,1,function(x) length(which(x==min(x,na.rm=T)))))-1
		mycut <- metabmiss * ncol(mymetab)
		keepers <- which(missnum < mycut)
		mymetab <- mymetab[keepers,]
		fmetab <- Biobase::AnnotatedDataFrame(data = Biobase::fData(inputData[["metabolite"]]))[keepers,]
	} else {
		print("No metabolite filtering by missing values is applied")
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

