#' Obtain lists of source IDs for inputting into RaMP or other pathway analysis
#' program
#' Requires KEGG ids and/or HMDB ids
#' Otherwise will output names given
#' @param inputResults results of IntLIM analysis
#' @param inputData IntLIM dataset containing feature data
#' @param outputMetab format of metabolites.  Either names or a list of all
#' metabolite ids.  'id' means source ID, 'name' means metabolite entry names,
#' 'mapping' provides names and source IDs
#' @return a list of source IDs or names for metabolites of interests
#'
#' @export
#'

getMetabRes <- function(inputResults, inputData, outputMetab='id'){

  fData.metab <- Biobase::fData(inputData[["metabolite"]])
  metab.list <- as.character(unique(inputResults@filt.results$metab))
  cols.fData.metab <- colnames(fData.metab)
  kegg.index <- grep('kegg',tolower(cols.fData.metab))
  hmdb.index <- grep('hmdb', tolower(cols.fData.metab))

  if((length(kegg.index) == 0 & length(hmdb.index)==0) | outputMetab == 'name' ){

    if(length(kegg.index) == 0 & length(hmdb.index)==0){
    print("No KEGG or HMDB columns in metabolite feature data.  Will print out metabolite names")
    }
    getMetabRes <- metab.list
    return(getMetabRes)
  }else{
  KEGG.id <- as.character(fData.metab[metab.list,kegg.index])
  HMDB.id <- as.character(fData.metab[metab.list,hmdb.index])
  RaMP.input <- KEGG.id
  RaMP.input[which(is.na(RaMP.input))] <- HMDB.id[which(is.na(RaMP.input))]

  metab.res <- data.frame('metab'=metab.list, KEGG.id, HMDB.id, RaMP.input)
  metab.id.list.obj <- as.character(metab.res$RaMP.input)
  metab.id.list.obj <- metab.id.list.obj[which(!is.na(metab.id.list.obj))]

  if (outputMetab == 'id'){
  getMetabRes <- metab.id.list.obj
  }else{
    getMetabRes <- metab.res
  }
  return(getMetabRes)
  }
}










#' Obtain lists of source IDs for inputting into RaMP or other pathway analysis
#' program
#' Requires sourceID for gene suitable for RaMP
#' Otherwise will output names given
#' @param inputResults results of IntLIM analysis
#' @param inputData IntLIM dataset containing feature data
#' @param outputGene format of metabolites.  Either names or a list of all
#' gene ids.  'id' means source ID, 'name' means gene entry names,
#' 'mapping' provides names and source IDs
#' @return a list of source IDs or names for genes of interests
#'
#' @export
#'

getGeneRes <- function(inputResults, inputData, outputGene='id'){

  fData.gene <- Biobase::fData(inputData[["expression"]])
  colnames.fData.gene <- colnames(fData.gene)
  sourceid.in <- 'sourceID' %in% colnames.fData.gene

  gene.list <- as.character(unique(inputResults@filt.results$gene))

    if(!sourceid.in | outputGene == 'name'){
    if(!sourceid.in){
      print('No column labeled sourceID in feature id.  Printing gene names')
    }
    getGeneRes <- gene.list
    return(getGeneRes)
  }else{

  ensembl.id.collec <- as.character(fData.gene[gene.list,'sourceID'])
  gene.id.list <- paste(ensembl.id.collec, collapse = ',')
  gene.id.list.obj <- unique(unlist(strsplit(gene.id.list, split = ',')))

  if(outputGene == 'id'){
  gene.id.list.obj <- gene.id.list.obj[which(!is.na(gene.id.list.obj))]
  getGeneRes <- gene.id.list.obj
  return(getGeneRes)
  }else{
  gene.res <- fData.gene[gene.list,c('id','sourceID')]
  getGeneRes <- gene.res
  return(getGeneRes)
  }

  }

}
