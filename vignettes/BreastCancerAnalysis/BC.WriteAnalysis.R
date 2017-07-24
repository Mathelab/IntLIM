library(readxl)
## Functions to assist

row.sds <- function(matrix.input){
    row.sds <- apply(matrix.input, 1, sd)
    names(row.sds) <- rownames(matrix.input)
    return(row.sds)
}

col.sds <- function(matrix.input){
    col.sds <- apply(matrix.input, 2, sd)
    names(col.sds) <- colnames(matrix.input)
    return(col.sds)
}

row.mean <- function(matrix.input){
    row.mean <- apply(matrix.input, 1, mean)
    names(row.mean) <- rownames(matrix.input)
    return(row.mean)
}

col.mean <- function(matrix.input){
    col.mean <- apply(matrix.input, 2, mean)
    names(col.mean) <- colnames(matrix.input)
    return(col.mean)
}

row.median <- function(matrix.input){
    row.median <- apply(matrix.input, 1, median)
    names(row.median) <- rownames(matrix.input)
    return(row.median)
}

col.median <- function(matrix.input){
    col.median <- apply(matrix.input, 2, mean)
    names(col.median) <- colnames(matrix.input)
    return(col.median)
}



## Import in gene expression data
## geneData is the original imported data file
geneData <- read.csv(file = 'geneCore.csv', row.names = 1)
## feature.annotation is the original annotation file
feature.annotation <- read.csv(file = 'annotation.csv')
## an X is attached to the probe ID
probe.id <- unlist(lapply(1:33297, function(x){return(paste('X', as.character(feature.annotation[x, 1]), sep = ""))}))
identical(probe.id, colnames(geneData))

##new fData.orig file - original fData object
fData.orig <- data.frame(probe.id, 'genesymbols' = feature.annotation$genesymbols)

## Remove gene probes that correspond to <NA>

## transpose of geneData data frame
geneData.t <- t(geneData)
no.na.vector <- which(is.na(feature.annotation$genesymbols) == FALSE)
geneData.t.noNA <- geneData.t[no.na.vector,]
fData.orig.noNA <- fData.orig[no.na.vector,]

## Identify duplicate probes and remove them

## List of unique.gene.symbols
unique.gene.symbols <- unique(fData.orig.noNA$genesymbols)

## new vector for number of each unique gene symbol
num.gene.symbols <- c()

# index of all genes to keep
keep.index <- c()

for (i in 1:length(unique.gene.symbols)){
    num.gene.symbols[i] <- length(which(fData.orig.noNA$genesymbols == unique.gene.symbols[i]))
    
    if (num.gene.symbols[i] == 1){
        keep.index <- c(keep.index, which(fData.orig.noNA$genesymbols == unique.gene.symbols[i]))
    }else{
        dup.indices <- which(fData.orig.noNA$genesymbols == unique.gene.symbols[i])
        
            dup.mat <- geneData.t.noNA[dup.indices,]
            dup.means <- row.mean(dup.mat)
            max.index <- which(dup.means == max(dup.means))[1]
            max.probe <- rownames(dup.mat)[max.index]
            keep.index <- c(keep.index, which(fData.orig.noNA$probe.id == max.probe))
        
    }
}

geneData.removeDup <- geneData.t.noNA[keep.index,]
fData.orig.removeDup <- fData.orig.noNA[keep.index,]
rownames(fData.orig.removeDup) <- c()

geneData.final <- geneData.removeDup
rownames(geneData.final) <- fData.orig.removeDup$genesymbols
zeros.vector <- rep(0, length(unique.gene.symbols))
fData.final <- data.frame('name' = fData.orig.removeDup$genesymbols, 
                          'id' = fData.orig.removeDup$genesymbols, 
                          'start' = zeros.vector, 'end' = zeros.vector, 
                          'chromosome' = zeros.vector)

write.csv(geneData.final, 'geneData.csv')
write.csv(fData.final, 'fData.gene.csv', row.names = FALSE)

## Extracting the Metabolon file

getMetabolon <- function(path, logmetab = FALSE, logbase = 2, phen.id.label = 'LHC', feature.id.label = 'BIOCHEMICAL'){
    #	path <- "metabolonfile.xls"
    
    #sampimp contains the excel read from the file path
    #Takes from the second sheet of the excel file (the imputed data)
    sampimp <- read_excel(path, sheet = 2, col_names = FALSE, col_types = NULL, na = "",skip = 0)
    sampimp <- as.data.frame(sampimp)
    #collect dimensions of the sampimp object
    dimSamp <- dim(sampimp)
    rows.samp <- dimSamp[1]
    cols.samp <- dimSamp[2]
    
    #look for the center cell with the HMDB_ID label
    #This will be pending discussion with a Metabolon representative
    #This is referred to as the 'index'
    get.index.col <- grep("HMDB_ID", sampimp)
    get.index.row <- grep("HMDB_ID", sampimp[,get.index.col])
    
    #index is actually HMDB_ID + another marker.  We should isolate that marker
    index.name.step1 <- unlist(strsplit(sampimp[get.index.row,get.index.col], "HMDB_ID"))
    index.name <- trimws(index.name.step1)
    
    ##get phenData
    #We will extract the phenoData from the file that we have extracted
    
    #pData contains the phenoData from the file
    pData <- sampimp[1:get.index.row, (get.index.col + 1):cols.samp]
    
    #the index column rows forms the row names for the pData
    rownames(pData) <- sampimp[1:get.index.row, get.index.col]
    #index.name is the name of the last row
    rownames(pData)[get.index.row] <- index.name
    
    #finds an entry in the pData rownames that contains an entry of 'id'
    #Will be used to set the id vector for the phenoData
    id.row.index <- which(rownames(pData) == phen.id.label)
    #print(phen.id.label)
    id <- pData[id.row.index,]
    #print(id)
    id <- as.integer(id)
    id <- as.character(id)
    
    #creates a pData data frame that contains pData + an id vector
    p.id <- id
    colnames(pData) <- id
    pData <- data.frame(t(pData))
    pData <- data.frame(id,pData)
    
    ##get featureData
    #We will extract the featureData from the file that we have extracted
    
    #fData contains the featureData from the file
    fData <- sampimp[(get.index.row+1):rows.samp,1:get.index.col]
    
    #The index row forms the column names for the fData
    col.fData <- dim(fData)[2]
    colnames(fData) <- sampimp[get.index.row, 1:get.index.col]
    #"HMDB_ID" is the name of the last column of featureData
    colnames(fData)[col.fData] <- "HMDB_ID"
    
    #METABOLON_ID is seen as the column for metabolite ids
    id.col.index <- which(colnames(fData) == feature.id.label)
    id <- fData[,id.col.index]
    #id <- as.integer(id)
    id <- as.character(id)
    f.id <- id
    rownames(fData) <- id
    
    #new data frame is created for the feature data
    fData <- data.frame(fData)
    fData <- data.frame(id,fData)
    
    #parse out the metabolite data matrix into metabData and convert into a matrix format
    metabData <- sampimp[(get.index.row+1):rows.samp, (get.index.col + 1):cols.samp]
    metabData <- as.matrix(metabData)
    orig.metabData <- metabData
    #Convert into a numeric matrix
    dim.md <- dim(metabData)
    metabData.number <- as.numeric(metabData)
    metabData <- matrix(metabData.number, dim.md[1], dim.md[2])
    
    #rownames and colnames of metabData are set as ids for phenoData and featureData
    colnames(metabData) <- p.id
    rownames(metabData) <- f.id
    
    
    if (logmetab == TRUE){
        metabData <- log(metabData,logbase)
    }
    
    
    getMetabolon <- list(metabData, pData, fData, orig.metabData)
    
    return(getMetabolon)
    
}

list.mb <- getMetabolon('metabolon.bc.data.xls')

metabData <- list.mb[[1]]
pData.metab <- list.mb[[2]]
fData.metab <- list.mb[[3]]

lhc.names <- rownames(pData.metab)

len_lhc <- length(lhc.names)

lhc.newnames <- c()

lhc.newnames <- unlist(lapply(lhc.names, function(x){return(paste('LHC', x, sep = ''))}))

colnames(metabData) <- lhc.newnames
rownames(pData.metab) <- lhc.newnames
pData.metab$LHC <- lhc.newnames
pData.metab$id <- lhc.newnames
metabData <- log2(metabData)
write.csv(pData.metab, "pData.csv", row.names = FALSE)
write.csv(metabData, "metabData.csv")
write.csv(fData.metab, "fData.metab.csv", row.names = FALSE)

