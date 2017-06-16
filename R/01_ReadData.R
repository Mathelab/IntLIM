#' Read in CSV file
#'
#' The metadata associated with data files to be analyzed in IntLim is supplied
#' as a CSV file with two columns and 6 rows: 
#'    type,filenames
#'    metabData,myfilename
#'    geneData,myfilename
#'    metabMetaData,myfilename (optional)
#'    geneMetaData,myfilename (optional)
#'    sampleMetaData,myfilename
#'
#' Note that all files supplied in the CSV file, and the CSV file itself should be placed in the same folder.  The software assumes will automatically retrieve the file path of
#' the input files (based on location of CSV files).  
#' Note also that the input data files should be in a specific format:
#'	metabData: rows are metabolites, columns are samples
#'	geneData: rows are genes, columns are samples
#'	metabMetaData: rows are metabolites, features are columns
#'	geneMetaData: rows are genes, features are columns
#'	sampleMetaData: rows are samples, features are columns
#' In addition, the first column of the sampleMetaData file is assumed to be the sample id, 
#' and those sample ids should match the columns of metabData and geneData (e.g. it is required
#' that all sample ids in the metabData and geneData are also in the sampleMetaDatafile).
#'
#'	
#' @param inputFile input file in CSV format (see Despcription)
#' @param metabid name of column from metabolite meta data to be used as id
#'      (required, must match metabolite abundances data)
#' @param geneid name of column from gene meta data to be used as id
#'	(required, must match gene expression data)
#' @param logmetab whether or not to log metabolite values (T/F)
#' @param loggene whether or not to log gene values (T/F)
#' @return MultiDataSet object with input data
#'
#' @examples
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "test.csv")
#' mydata <- ReadData(csvfile,metabid='BIOCHEMICAL',geneid='X')
#' @export
ReadData <- function(inputFile,metabid,geneid, logmetab=FALSE,loggene=FALSE){
    # Check that file exists
    if (!file.exists(inputFile)) {
        stop("CSV input file does not exist")
    }
    # Make into df to make access easier
    csvfile <- as.data.frame(utils::read.csv(inputFile, header=TRUE,row.names=1))

    # Check column names are correct
    if (colnames(csvfile)!="filenames") {
        stop("Check column names of input files.  'type' and 'filenames' are required")
    }
   
    # Check that all types required are present
    mytypes <- c("metabData","geneData","metabMetaData","geneMetaData",
         "sampleMetaData")
    mymatches <- as.numeric(lapply(mytypes,function(x) 
		length(which(rownames(csvfile)==x))))
    if(sum(mymatches)!=5) {
	stop("The column 'type' contains non-allowed entries (See Description)")}
 
    # Check that files exist then read them in one by one
    if(!file.exists(as.character(csvfile['metabData',]))) {
	stop(paste("File", csvfile['metabData',], "does not exist"))} else {
    MData<-utils::read.csv(as.character(csvfile['metabData',]),row.names = 1)}

    if(!file.exists(as.character(csvfile['geneData',]))) {
        stop(paste("File", csvfile['geneData',], "does not exist"))} else {
    GData<-utils::read.csv(as.character(csvfile['geneData',]),row.names = 1)}

    if(as.character(csvfile['metabMetaData',])=="") {
	warning("No metadata provided for metabolites");MmetaData<-NULL; } else if
    (!file.exists(as.character(csvfile['metabMetaData',]))) {
        stop(paste("File", csvfile['metabMetaData',], "does not exist"))} else {
    MmetaData<-utils::read.csv(as.character(csvfile['metabMetaData',]))}

    if(as.character(csvfile['geneMetaData',])=="") {
        warning("No metadta provided for genes");GmetaData<-NULL;} else if
    (!file.exists(as.character(csvfile['geneMetaData',]))) {
        stop(paste("File", csvfile['geneMetaData',], "does not exist"))} else {
    GmetaData<-utils::read.csv(as.character(csvfile['geneMetaData',]))}

    if(!file.exists(as.character(csvfile['sampleMetaData',]))) {
        stop(paste("File", csvfile['sampleMetaData',], "does not exist"))} else {
    pData<-utils::read.csv(as.character(csvfile['sampleMetaData',]),row.names = 1)}

    #Create Multi
    GMdata <- CreateIntLimObject(genefdata=GmetaData, metabfdata=MmetaData,
	metabid=metabid, geneid=geneid,pdata=pData,metabdata=MData,
	genedata=GData,logmetab=logmetab,loggene=loggene)
	
	#genefdata=GmetaData; metabfdata=MmetaData;metabid="BIOCHEMICAL";geneid="X";
	#pdata=pData;metabdata=MData;genedata=GData;logmetab=F;loggene=F

    print("CreateMultiDataSet created")
    return(GMdata)
}



