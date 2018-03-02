#' Generic function to create constructor of MultiDataSet gene object
#'
#' @include MetaboliteSet_addMetabolite.R
#' @include AllClasses.R 
#'
#' @param genefdata gene meta data
#' @param metabfdata metabolite meta data
#' @param pdata sample meta data
#' @param geneid name of column from metabolite meta data to be used as id 
#'	(required if a gene meta data file is present, must match gene expression matrix))
#' @param metabid name of column from gene meta data to be used as id
#'      (required if a metabolite meta data file is present, must match metabolite abundances matrix))
#' @param metabdata metabolite abundances (samples are in columns)
#' @param genedata gene expression (samples are in columns)
#' @param logmetab T/F 
#' @param loggene T/F
CreateIntLimObject <- function(genefdata, metabfdata, pdata, geneid, metabid, 
	metabdata, genedata, logmetab=FALSE,loggene=FALSE) {

	# Check that feature data and abundance data metabolites corresponds
        if (!is.null(metabfdata)) {
        if(length(which(colnames(metabfdata)=='id'))!=1) {
                stop(paste("metabid provided",metabid,"does not exist in metabolite meta data file"))} else if 
	(length(intersect(rownames(metabdata),as.character(metabfdata[,metabid])))<nrow(metabdata)){
                stop("Metabolites in abundance data file and metabolite meta data file are not equal")} else {
                myind <- as.numeric(lapply(rownames(metabdata),function(x) {
                        which(as.character(metabfdata[,'id'])==x)[1]}))
                        metabpdata<-pdata[myind,]}

	rownames(metabfdata)=as.character(metabfdata[,'id'])
        }

        # Check that samples data and abundance data samples correspond
        if(length(intersect(colnames(metabdata),rownames(pdata)))<ncol(metabdata)){
                stop("All samples in abundance data file must be in metabolite meta data file")
        } else {
                myind <- as.numeric(lapply(colnames(metabdata),function(x) {
                        which(rownames(pdata)==x)[1]}))
                        metabpdata<-pdata[myind,]
        }

	#new data frames are set for phenoData and featureData
	metabpdata$id=rownames(metabpdata)
	metabphenoData <- Biobase::AnnotatedDataFrame(data = metabpdata)

	if (logmetab == TRUE){
		metabdata <- log2(metabdata)
	}
	
	if(is.null(metabfdata)) {
		metabfdata <- data.frame(id = rownames(metabdata),stringsAsFactors=FALSE)
                rownames(metabfdata) <- metabfdata[,1]
	}
	# Make sure order of feature data is the same as the data matrix:
	myind=as.numeric(lapply(rownames(metabdata),function(x) which(metabfdata[,"id"]==x)))
        metabfdata <- data.frame(metabfdata[myind,],stringsAsFactors=FALSE)
	rownames(metabfdata) <- metabfdata[,1]
	metabfeatureData <- Biobase::AnnotatedDataFrame(data = metabfdata)
	metab.set <- methods::new("MetaboliteSet",metabData = metabdata, 
		phenoData = metabphenoData, featureData = metabfeatureData)
	

	#####  Now the genes
	# Check that feature data and gene expression data corresponds
       if(!is.null(genefdata)) {
       if(length(which(colnames(genefdata)=="id"))!=1) {
                stop(paste("geneid provided",geneid,"does not exist in gene meta data file"))
        } else if(length(intersect(rownames(genedata),as.character(genefdata[,'id']))) < nrow(genedata)){
                stop("Genes in expression data file and gene meta data file are not equal")
        } else {
                myind <- as.numeric(lapply(rownames(genedata),function(x) {
                        which(as.character(genefdata[,'id'])==x)[1]}))
                        genepdata<-pdata[myind,]
        }

        rownames(genefdata)=as.character(genefdata[,'id'])
        }
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
                genedata <- log2(genedata)
        }

        gene.set <- Biobase::ExpressionSet(assayData=as.matrix(genedata))
	if(is.null(genefdata)) {
		genefdata <- data.frame(id = rownames(genedata))
		rownames(genefdata) <- genefdata[,1]
        }
	if(length(which(colnames(genefdata)=="chromosome"))==0) {
		genefdata$chromosome <- rep("chr",nrow(genefdata))}
	if(length(which(colnames(genefdata)=="start"))==0) {
	        genefdata$start <- rep(0,nrow(genefdata))}
	if(length(which(colnames(genefdata)=="end"))==0) {
	        genefdata$end <- rep(0,nrow(genefdata))}
	# Make sure that the order of genefdata is the same as the input data
	myind=as.numeric(lapply(rownames(genedata),function(x) which(genefdata$id==x)))
	genefdata <- genefdata[myind,]
	Biobase::fData(gene.set) <- data.frame(genefdata,stringAsFactors=FALSE)
	genepdata$id=rownames(genepdata)
	Biobase::pData(gene.set) <- genepdata
	

	multi <- MultiDataSet::createMultiDataSet()
	multi1 <- MultiDataSet::add_genexp(multi, gene.set)
	multi2 <- add_metabolite(multi1, metab.set)
	return(multi2)
}


#' Function that returns a list of all data for samples in common between metabolite and gene datasets
#'
##' @include AllClasses.R
#'
#' @import MultiDataSet
#' @param inputData MultiDataSet object (output of ReadData())
#' @param stype category to color-code by (can be more than two categories)
#' @param covar vector of additional variables to be incorporated into model
#' @param class.covar class of additional variables
getCommon <- function(inputData,stype=NULL, covar = NULL, class.covar = NULL) {
   incommon<-MultiDataSet::commonSamples(inputData)
   mp <- Biobase::pData(incommon[["metabolite"]])
   gp <- Biobase::pData(incommon[["expression"]])

   if(all.equal(mp[,stype],gp[,stype])[1] != TRUE) {
        stop(paste("The column", stype,"for the samples in common between the metabolite and gene datasets are not equal.  Please check your input."))
    }


   gene <- Biobase::assayDataElement(inputData[["expression"]], 'exprs')
   metab <- Biobase::assayDataElement(inputData[["metabolite"]], 'metabData')

   # Force the order to be the same, in case it isn't
   p <- mp[rownames(gp),]
   p.0 <- p
   metab <- metab[,colnames(gene)]

   if(!is.null(stype)) {
	p <- p[,stype]
        uniqp <- unique(p)

	uniqtypes <- unique(p)
        # Deal with missing values or ""
        if(length(which(p==""))>0) {
                new.p <- p[which(p!="")]
		metab <- metab[,which(p!="")]
		gene <- gene[,which(p!="")]
		p <- new.p
        }
        if(length(which(is.na(p)))>0) {
                new.p <- p[which(!is.na(p))]
                metab <- metab[,which(!is.na(p))]
                gene <- gene[,which(!is.na(p))]
		p <- new.p
        }
   }
   
   if(!is.null(covar)){
       
       if (length(covar %in% colnames(p.0)) != sum(covar %in% colnames(p.0))){
           stop("Additional variable names not in pData")
       }
       covar_matrix <- p.0[colnames(gene),covar, drop = FALSE]
       na.covar <- which(is.na(covar_matrix) | covar_matrix == '',arr.ind = TRUE)
       na.covar.list <- unique(rownames(na.covar))
       new.overall.list <- setdiff(colnames(gene), na.covar.list)
       
       covar_matrix <- covar_matrix[new.overall.list,,drop = FALSE]
       
       class.var <- apply(covar_matrix,2,class)
       
       gene <- gene[,new.overall.list]
       metab <- metab[,new.overall.list]
       p <- p.0[new.overall.list,stype]
       
       if(!(is.null(class.covar))){
           
           if(length(class.covar) != length(covar)){
               stop("lengths of covar and class.covar not the same")
           }
           len.covar <- length(covar)
           for(i in 1:len.covar){
               if(class.covar[i] == 'numeric'){
                   
                   covar_matrix[,i] <- as.numeric(covar_matrix[,i])
                   
               }else{
                   
                   covar_matrix[,i] <- as.factor(as.character(covar_matrix[,i]))
                   
               }
           }
       }
       
   }else{
       covar_matrix <- NULL
   }

   # Check that everything is in right order
   if(!all.equal(rownames(mp),rownames(gp)) || !all.equal(colnames(metab),colnames(gene))){ 
	stop("Something went wrong with the merging!  Sample names of input files may not match.")
   } else {
   out <- list(p=as.factor(as.character(p)),gene=gene,metab=metab, covar_matrix=covar_matrix)
   }
   return(out)
}

#' Function that runs linear models and returns interaction p-values.
#'
#' @include MetaboliteSet_addMetabolite.R
#' @include AllClasses.R
#'
#' @param incommon MultiDataSet object (output of ReadData()) with gene
#' @param outcome 'metabolite' or 'gene' must be set as outcome/independent variable 
#' (default is 'metabolite')
#' @param type vector of sample type (by default, it will be used in the interaction term).
#' Only 2 categories are currently supported.
#' @param covar vector of additional vectors to consider
RunLM <- function(incommon, outcome="metabolite", type=NULL, covar=NULL) { 

    gene <- incommon$gene
    metab <- incommon$metab
 
    uniqtypes <- unique(type)
    if(length(uniqtypes)!=2) {
	stop("The number of unique categores is not 2.")
    }

    genesd1 <- as.numeric(apply(gene[,which(type==uniqtypes[1])],1,function(x) stats::sd(x,na.rm=T)))
    metabsd1 <- as.numeric(apply(metab[,which(type==uniqtypes[1])],1,function(x) stats::sd(x,na.rm=T)))
    genesd2 <- as.numeric(apply(gene[,which(type==uniqtypes[2])],1,function(x) stats::sd(x,na.rm=T)))
    metabsd2 <- as.numeric(apply(metab[,which(type==uniqtypes[2])],1,function(x) stats::sd(x,na.rm=T)))

    mymessage=""
    if(length(which(genesd1==0))>0 || length(which(genesd2==0))>0) {
	toremove <- c(which(genesd1==0),which(genesd2==0))
	gene <- gene[-toremove,]
	mymessage <- c(mymessage,paste("Removed",length(toremove),"genes that had a standard deviation of 0:"))
	mymessage <- c(mymessage,rownames(gene)[toremove])
    }
    if(length(which(metabsd1==0))>0 || length(which(metabsd2==0))>0) {
        toremove <- c(which(metabsd1==0),which(metabsd2==0))
        metab <- metab[-toremove,]
        mymessage <- c(mymessage,paste("Removed",length(toremove),"metabolites that had a standard deviation of 0:"))
        mymessage <- c(mymessage,rownames(metab)[toremove])
    }

    if (outcome == "metabolite") {
        arraydata <- data.frame(metab)
        #form <- stats::formula(m ~ g + type + g:type)
        # Retrieve pvalues by iterating through each gene
        numgenes <- nrow(gene)
	numprog <- round(numgenes*0.1)
	form.add <- "Y ~ g + type + g:type"
	    if (!(is.null(covar))){
	        form.add <- "Y ~ g + type + g:type"
	        
	        len.covar <- length(covar)
	        for (i in 1:len.covar){
	            form.add <- paste(form.add, '+', covar[i])
	        }
	    }
	
	
        list.pvals <- lapply(1:numgenes, function(x) {
                #print(x)
                g <- as.numeric(gene[x,])
                
                if(is.null(covar)){
                clindata <- data.frame(g, type)
                }else{
                    clindata <- data.frame(g, type, incommon$covar_matrix)
                }
                
                mlin <- getstatsOneLM(stats::as.formula(form.add), clindata = clindata,
                        arraydata = arraydata)
                term.pvals <- rownames(mlin$p.value.coeff)
                index.interac <- grep('g:type', term.pvals)
                
                p.val.vector <- as.vector(mlin$p.value.coeff[index.interac,])
                
                
                #p.val.vector <- as.vector(mlin@p.value.coeff[4,])
                # Print out progress every 1000 genes
                if (x %% numprog == 0){
                    progX <- round(x/numgenes*100)
                    print(paste(progX,"% complete"))
                }
                return(p.val.vector)
        })
    mat.pvals <- do.call(rbind, list.pvals)
    # adjust p-values
    row.pvt <- dim(mat.pvals)[1]
    col.pvt <- dim(mat.pvals)[2]
    myps <- as.vector(mat.pvals)
    mypsadj <- stats::p.adjust(myps, method = 'fdr')
    mat.pvalsadj <- matrix(mypsadj, row.pvt, col.pvt)

    rownames(mat.pvals) <- rownames(mat.pvalsadj) <- rownames(gene)
    colnames(mat.pvals) <- colnames(mat.pvalsadj) <- rownames(metab)
    } else if (outcome == "gene") {
        arraydata <- data.frame(gene)
        #form <- stats::formula(g ~ m + type + m:type)

        # Retrieve pvalues by iterating through each gene
        nummetab <- nrow(metab)
	numprog <- round(nummetab*0.1)
	
	form.add <- "Y ~ m + type + m:type"
	if (!(is.null(covar))){
	    form.add <- "Y ~ m + type + m:type"
	    
	    len.covar <- length(covar)
	    for (i in 1:len.covar){
	        form.add <- paste(form.add, '+', covar[i])
	    }
	}
        list.pvals <- lapply(1:nummetab, function(x) {
                #print(x)
                m <- as.numeric(metab[x,])
                
                if(is.null(covar)){
                    clindata <- data.frame(m, type)
                }else{
                    clindata <- data.frame(m, type, incommon$covar_matrix)
                }
               
                mlin <- getstatsOneLM(stats::as.formula(form.add), clindata = clindata,
                        arraydata = arraydata)
                
                term.pvals <- rownames(mlin$p.value.coeff)
                index.interac <- grep('m:type', term.pvals)
                p.val.vector <- as.vector(mlin$p.value.coeff[index.interac,])
                #p.val.vector <- as.vector(mlin@p.value.coeff[4,])
                # Print out progress every 1000 genes
                if (x %% numprog == 0){
                    progX <- round(x/nummetab*100)
                    print(paste(progX,"% complete"))
                }
                return(p.val.vector)
        })
    mat.pvals <- do.call(rbind, list.pvals)
    # adjust p-values
    row.pvt <- dim(mat.pvals)[1]
    col.pvt <- dim(mat.pvals)[2]
    myps <- as.vector(mat.pvals)
    mypsadj <- stats::p.adjust(myps, method = 'fdr')
    mat.pvalsadj <- matrix(mypsadj, row.pvt, col.pvt)

    rownames(mat.pvals) <- rownames(mat.pvalsadj) <- rownames(metab)
    colnames(mat.pvals) <- colnames(mat.pvalsadj) <- rownames(gene)

    mat.pvals <- t(mat.pvals)
    } else {
        stop("outcome must be either 'metabolite' or 'gene'")
    }

    myres <- methods::new('IntLimResults', interaction.pvalues=mat.pvals,
		interaction.adj.pvalues = mat.pvalsadj,
		warnings=mymessage)
    return(myres)
}

#' Function that runs linear models for one gene vs all metabolites
#'
#' @include AllClasses.R
#'
#' @param form LM formulat (typically m~g+t+g:t)
#' @param clindata data frame with 1st column: expression of one analyte; 2nd column
#' sample type (e.g. cancer/non-cancer)
#' @param arraydata matrix of metabolite values
    getstatsOneLM <- function(form, clindata, arraydata) {
	call=match.call()
        YY <- t(arraydata)                      # the data matrix
        EY <- apply(YY, 2, mean)                # its mean vector
        SYY <- apply(YY, 2, function(y) {sum(y^2)}) - nrow(YY)*EY^2     # sum of squares after centering
        clindata <- data.frame(y=YY[,1], clindata)
        dimnames(clindata)[[2]][1] <- 'Y'
        X <- stats::model.matrix(form, clindata)       # contrasts matrix
        N = dim(X)[1]
        p <- dim(X)[2]
        XtX <- t(X) %*% X
        ixtx <- solve(XtX)
        bhat <- ixtx %*% t(X) %*% YY            # Use the pseudo-inverse to estimate the parameters
        yhat <- X %*% bhat                      # Figure out what is predicted by the model
        # Now we partition the sum-of-square errors
        rdf <- ncol(X)-1                        # number of parameters in the model
        edf <- nrow(YY)-rdf-1                   # additional degrees of freedom
        errors <- YY - yhat                     # difference between observed and model predictions
        sse <- apply(errors^2, 2, sum)  # sum of squared errors over the samples
        mse <- sse/edf                  # mean squared error
        ssr <- SYY - sse                        # regression error
        msr <- ssr/rdf                  # mean regression error
        fval <- msr/mse                 # f-test for the overall regression
        pfval <- 1-stats::pf(fval, rdf, edf)           # f-test p-values

        stderror.coeff <- sapply(mse,function(x){sqrt(diag(ixtx)*x)})
        t.coeff <- bhat/stderror.coeff
        p.val.coeff <- 2*stats::pt(-abs(t.coeff),df = (N-p))
        #methods::new('IntLimModel', call=call, model=form,
         list(# call=call, model=form,
                #coefficients=bhat,
                # predictions=yhat,
                #df=c(rdf, edf),
                #sse=sse,
                #ssr=ssr,
                #F.statistics=fval,
                #F.p.values=pfval
                #std.error.coeff = stderror.coeff,
                #t.value.coeff = t.coeff,
                p.value.coeff = p.val.coeff # interaction p-value
         )
        }


