#' Generic function to create constructor of MultiDataSet gene object
#'
#' @include MetaboliteSet_addMetabolite.R
#' @include AllClasses.R 
#'
#' @param genefdata gene meta data
#' @param metabfdata metabolite meta data
#' @param pdata sample meta data
#' @param geneid name of column from metabolite meta data to be used as id 
#'	(required, must match metabolite abundances matrix))
#' @param metabid name of column from gene meta data to be used as id
#'      (required, must match gene abundances matrix))
#' @param metabdata metabolite abundances (samples are in columns)
#' @param genedata gene expression (samples are in columns)
#' @param logmetab T/F 
#' @param loggene T/F

CreateIntLimObject <- function(genefdata, metabfdata, pdata, geneid, metabid, 
	metabdata, genedata, logmetab=FALSE,loggene=FALSE) {

	# Check that feature data and abundance data metabolites corresponds
        if(length(which(colnames(metabfdata)==metabid))!=1) {
                stop(paste("metabid provided",metabid,"does not exist in metabolite meta data file"))} else if 
	(length(intersect(rownames(metabdata),as.character(metabfdata[,metabid])))<nrow(metabdata)){
                stop("Metabolites in abundance data file and metabolite meta data file are not equal")} else {
                myind <- as.numeric(lapply(rownames(metabdata),function(x) {
                        which(as.character(metabfdata[,metabid])==x)[1]}))
                        metabpdata<-pdata[myind,]}

	rownames(metabfdata)=as.character(metabfdata[,metabid])

        # Check that samples data and abundance data samples correspond
        if(length(intersect(colnames(metabdata),rownames(pdata)))<ncol(metabdata)){
                stop("Samples in metabolite abundance data file and sample meta data file are not equal")
        } else {
                myind <- as.numeric(lapply(colnames(metabdata),function(x) {
                        which(rownames(pdata)==x)[1]}))
                        metabpdata<-pdata[myind,]
        }

	#new data frames are set for phenoData and featureData
	metabpdata$id=rownames(metabpdata)
	metabphenoData <- Biobase::AnnotatedDataFrame(data = metabpdata)
	metabfeatureData <- Biobase::AnnotatedDataFrame(data = metabfdata)

	if (logmetab == TRUE){
		metabdata <- log2(metabdata)
	}

	metab.set <- methods::new("MetaboliteSet",metabData = metabdata, 
		phenoData = metabphenoData, featureData = metabfeatureData)

	#####  Now the genes
	# Check that feature data and gene expression data corresponds
       if(length(which(colnames(genefdata)==geneid))!=1) {
                stop(paste("geneid provided",geneid,"does not exist in gene meta data file"))
        } else if(length(intersect(rownames(genedata),as.character(genefdata[,geneid]))) < nrow(genedata)){
                stop("Genes in expression data file and gene meta data file are not equal")
        } else {
                myind <- as.numeric(lapply(rownames(genedata),function(x) {
                        which(as.character(genefdata[,geneid])==x)[1]}))
                        genepdata<-pdata[myind,]
        }

        rownames(genefdata)=as.character(genefdata[,geneid])

        # Check that samples data and abundance data samples correspond
        if(length(intersect(colnames(genedata),rownames(pdata)))<ncol(genedata)){ 
                stop("Samples in expression data file and sample meta data file are not equal")
        } else {
		myind <- as.numeric(lapply(colnames(genedata),function(x) {
			which(rownames(pdata)==x)[1]}))
			genepdata<-pdata[myind,]
	}

        #new data frames are set for phenoData and featureData
        if (loggene == TRUE){
                metabdata <- log2(genedata)
        }

        gene.set <- Biobase::ExpressionSet(assayData=as.matrix(genedata))
	if(length(which(colnames(genefdata)=="chromosome"))==0) {
		genefdata$chromosome <- rep("chr",nrow(genefdata))}
	if(length(which(colnames(genefdata)=="start"))==0) {
	        genefdata$start <- rep(0,nrow(genefdata))}
	if(length(which(colnames(genefdata)=="end"))==0) {
	        genefdata$end <- rep(0,nrow(genefdata))}
	Biobase::fData(gene.set) <- genefdata
	genepdata$id=rownames(genepdata)
	Biobase::pData(gene.set) <- genepdata
	

	multi <- MultiDataSet::createMultiDataSet()
	multi1 <- MultiDataSet::add_genexp(multi, gene.set)

#    methods::setGeneric("add_metabolite", function(object, metabSet, warnings = TRUE, ...)
#    base::standardGeneric("add_metabolite")
#    )

 #   methods::setMethod(
 #       f = "add_metabolite",
 #       signature = c("MultiDataSet", "MetaboliteSet"),
 #       definition = function(object, metabSet, warnings = TRUE, ...) {
 #           ## Add given MetaboliteSet as 'metabolite'
 #           object <- MultiDataSet::add_eset(object, metabSet, dataset.type = "metabolite", GRanges = NA, ...)
 #           return(object)
 #       })

	multi2 <- add_metabolite(multi1, metab.set)
	return(multi2)
}




