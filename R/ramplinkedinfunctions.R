#' Obtain lists of source IDs for inputting into RaMP or other pathway analysis
#' program
#' Requires sourceIDs and IDtypes
#' Otherwise will output names given
#' The IDtype should be as below  \cr
#' +------------+ \cr
#'| IDtype     | \cr
#'  +------------+ \cr
#'  | CAS        | \cr
#'  | chebi      | \cr
#'  | chemspider | \cr
#'  | hmdb       | \cr
#'  | kegg       | \cr
#'  | LIPIDMAPS  | \cr
#'  | pubchem    | \cr
#'  +------------+ \cr
#' @param inputResults results of IntLIM analysis
#' @param inputData IntLIM dataset containing feature data
#' @param outputMetab format of metabolites.  Either names or a list of all
#' metabolite ids.  'id' means source ID, 'name' means metabolite entry names,
#' 'mapping' provides names and source IDs
#' @return a list of source IDs or names for metabolites of interests
#'
#' @export
#'

getMetabList <- function(inputResults, inputData, outputMetab='id'){

  fData.metab <- Biobase::fData(inputData[["metabolite"]])
  metab.list <- as.character(unique(inputResults@filt.results$metab))
  cols.fData.metab <- colnames(fData.metab)

  if (outputMetab == 'id' | outputMetab == 'mapping'){


  if('sourceID' %in% cols.fData.metab & 'IDtype' %in% cols.fData.metab){
    len.metabs <- nrow(fData.metab)
    mappinglist <- c()

    mappingIDs <- function(index){
      id.source <- as.character(fData.metab$IDtype[index])
      id.name <- as.character(fData.metab$sourceID[index])
      id.list <- unlist(strsplit(id.name,split=','))
      id.DB.added <- unlist(lapply(id.list, function(x){return(paste(id.source, ":", x, sep = ""))}))


      mappingIDs <- paste(id.DB.added, collapse = ',')

      if(is.na(id.source) | is.na(id.name)){
        mappingIDs <- NA
      }

      return(mappingIDs)
    }

    mappinglist <- unlist(lapply(1:nrow(fData.metab), mappingIDs))
    mapping.complete <- data.frame('name' = fData.metab$id, 'mapping' = mappinglist)
    rownames(mapping.complete) <- mapping.complete$name
    mapping.summary <- mapping.complete[metab.list,]

    if (outputMetab == 'id'){
      mapping.res.list.string <- paste(as.character(mapping.summary$mapping[!is.na(mapping.summary$mapping)],
                                                    collapse = ','))
      mapping.res.list <- unique(unlist(strsplit(mapping.res.list.string, split = ',')))
      getMetabList <- mapping.res.list[!is.na(mapping.res.list)]
      return(getMetabList)
    }else{ #asking for mapping
      getMetabList <- mapping.summary
      return(getMetabList)
    }
  }else{
    print("No sourceID or IDtype column.  Outputting only names")
    getMetabList <- metab.list
    return(getMetabList)
  }

  }else{
    getMetabList >- metab.list
    return(getMetabList)

  }
}



#' Obtain lists of source IDs for inputting into RaMP or other pathway analysis
#' program
#' Requires sourceIDs and IDtypes
#' Otherwise will output names given
#' The IDtype should be as below \cr
#' +--------------------+ \cr
#' | IDtype             | \cr
#'  +--------------------+ \cr
#'  | enzymeNomenclature | \cr
#'  | ensembl            | \cr
#'  | entrez             | \cr
#'  | hmdb               | \cr
#'  | kegg               | \cr
#'  | uniprot            | \cr
#'  +--------------------+ \cr
#' @param inputResults results of IntLIM analysis
#' @param inputData IntLIM dataset containing feature data
#' @param outputGene format of genes.  Either names or a list of all
#' gene ids.  'id' means source ID, 'name' means gene entry names,
#' 'mapping' provides names and source IDs
#' @return a list of source IDs or names for genes of interests
#'
#' @export
#'

getGeneList <- function(inputResults, inputData, outputGene='id'){

  fData.gene <- Biobase::fData(inputData[["expression"]])
  gene.list <- as.character(unique(inputResults@filt.results$gene))
  cols.fData.gene <- colnames(fData.gene)

  if (outputGene == 'id' | outputGene == 'mapping'){


    if('sourceID' %in% cols.fData.gene & 'IDtype' %in% cols.fData.gene){
      len.metabs <- nrow(fData.gene)
      mappinglist <- c()

      mappingIDs <- function(index){
        id.source <- as.character(fData.gene$IDtype[index])
        id.name <- as.character(fData.gene$sourceID[index])
        id.list <- unlist(strsplit(id.name,split=','))
        id.DB.added <- unlist(lapply(id.list, function(x){return(paste(id.source, ":", x, sep = ""))}))


        mappingIDs <- paste(id.DB.added, collapse = ',')

        if(is.na(id.source) | is.na(id.name)){
          mappingIDs <- NA
        }

        return(mappingIDs)
      }

      mappinglist <- unlist(lapply(1:nrow(fData.gene), mappingIDs))
      mapping.complete <- data.frame('name' = fData.gene$id, 'mapping' = mappinglist)
      rownames(mapping.complete) <- mapping.complete$name
      mapping.summary <- mapping.complete[gene.list,]

      if (outputGene == 'id'){
        mapping.res.list.string <- paste(as.character(mapping.summary$mapping[!is.na(mapping.summary$mapping)],
                                                      collapse = ','))
        mapping.res.list <- unique(unlist(strsplit(mapping.res.list.string, split = ',')))
        getGeneList <- mapping.res.list[!is.na(mapping.res.list)]
        return(getGeneList)
      }else{ #asking for mapping
        getGeneList <- mapping.summary
        return(getGeneList)
      }
    }else{
      print("No sourceID or IDtype column.  Outputting only names")
      getGeneList <- gene.list
      return(getGeneList)
    }

  }else{
    getGeneList >- gene.list
    return(getGeneList)

  }
}


