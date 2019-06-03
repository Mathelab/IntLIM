#' Retrieve significant gene-metabolite pairs, based on adjusted p-values.
#' For each gene-metabolite pair that is statistically significant, calculate the
#' correlation within group1 (e.g. cancer) and the correlation within group2 (e.g.
#' non-cancer).  Users can then remove pairs with a difference in correlations between
#' groups 1 and 2 less than a user-defined threshold.
#'
#' @include internalfunctions.R
#'
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param inputData MultiDataSet object (output of ReadData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param pvalcutoff cutoff of FDR-adjusted p-value for filtering (default 0.05)
#' @param diffcorr cutoff of differences in correlations for filtering (default 0.5)
#' @param corrtype spearman or pearson or other parameters allowed by cor() function (default
#' spearman)
#' @param treecuts user-selected number of clusters (of gene-metabolite pairs) to cut the tree into
#' @return IntResults object with model results (now includes correlations)
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata,treecuts=2)
#' }
#' @export
ProcessResults <- function(inputResults,
				inputData,
				pvalcutoff=0.05,
				diffcorr=0.5,
				corrtype="spearman",
				treecuts = 0
			){

	if(inputResults@outcome == "metabolite") {
		mydat <-inputResults@interaction.adj.pvalues}
		#mydat <- reshape2::melt(inputResults@interaction.adj.pvalues)}
	else if (inputResults@outcome == "gene") {
		mydat <-t(inputResults@interaction.adj.pvalues)}
                #mydat <- reshape2::melt(t(inputResults@interaction.adj.pvalues))}

	incommon <- getCommon(inputData,inputResults@stype)
	p <- incommon$p
	gene <- incommon$gene
	metab <- incommon$metab
	if(length(unique(p)) !=2) {
 		stop(paste("IntLim currently requires only two categories.  Make sure the column",inputResults@stype,"only has two unique values"))
    }

	gp1 <- which(p == unique(p)[1])
	cor1.m <- stats::cor(t(gene[rownames(mydat),gp1]),t(metab[colnames(mydat),gp1]),method=corrtype)
        gp2 <- which(p == unique(p)[2])
        cor2.m <- stats::cor(t(gene[rownames(mydat),gp2]),t(metab[colnames(mydat),gp2]),method=corrtype)

	if(pvalcutoff == 1) { #(no filtering)
		temp <- reshape2::melt(cor1.m)
		fincor1 <- as.numeric(temp[,"value"])
		temp <- reshape2::melt(cor2.m)
		fincor2 <- as.numeric(temp[,"value"])
		genenames <- as.character(temp[,1])
		metabnames <- as.character(temp[,2])
	} else {
		keepers <- which(mydat <= pvalcutoff, arr.ind=T)
		fincor1 <- as.numeric(apply(keepers,1,function(x)
			cor1.m[x[1],x[2]]))
                fincor2 <- as.numeric(apply(keepers,1,function(x)
                        cor2.m[x[1],x[2]]))
		genenames <- as.character(rownames(cor1.m)[keepers[,1]])
		metabnames <- as.character(colnames(cor1.m)[keepers[,2]])
	}

        mydiffcor = abs(fincor1-fincor2)

	keepers2 <- which(mydiffcor >= diffcorr)

	inputResults@filt.results <- data.frame(metab=metabnames[keepers2],
		gene=genenames[keepers2])
	inputResults@filt.results <- cbind(inputResults@filt.results,fincor1[keepers2],fincor2[keepers2])
	colnames(inputResults@filt.results)[3:4]=paste0(setdiff(as.character(unlist(unique(p))),""),"_cor")

	diff.corr <- inputResults@filt.results[,4] - inputResults@filt.results[,3]

	inputResults@filt.results <- cbind(inputResults@filt.results, diff.corr)
	if(inputResults@outcome == "metabolite") {
                adjp <- reshape2::melt(inputResults@interaction.adj.pvalues)
		p <-  reshape2::melt(inputResults@interaction.pvalues)
	} else if (inputResults@outcome == "gene") {
                adjp <- reshape2::melt(t(inputResults@interaction.adj.pvalues))
		p <- reshape2::melt(t(inputResults@interaction.pvalues))
	}

	cornames <- paste(as.character(inputResults@filt.results[,"metab"]),as.character(inputResults@filt.results[,"gene"]))
	rownames(p) <- paste(as.character(p[,2]),as.character(p[,1]))
	rownames(adjp) <- paste(as.character(adjp[,2]),as.character(adjp[,1]))
	outp <- p[cornames,]
	outpadj <- adjp[cornames,]

	inputResults@filt.results = cbind(inputResults@filt.results,outp$value, outpadj$value)
	colnames(inputResults@filt.results)[6:7]=c("Pval","FDRadjPval")


	if (treecuts > 0){

	hc.rows<- stats::hclust(stats::dist(inputResults@filt.results[,c(3,4)]))
	cluster <- stats::cutree(hc.rows, k = treecuts)

	inputResults@filt.results = cbind(inputResults@filt.results, cluster)


	}

print(paste(nrow(inputResults@filt.results), 'gene-metabolite pairs found given cutoffs'))
return(inputResults)
}

#' Retrieve significant gene-metabolite pairs (aka filter out nonsignificant pairs) based on value of gene:type interaction coefficient from linear model
#' @param inputResults IntLimResults object with model results (output of RunIntLim())
#' @param InteractionCoeffcutoff Smallest interaction coefficient that will be graphed (positive or negative)
#' @return IntLimResults object with model results (now includes filt.results data)
#' @export
ProcessResultsContinuous<- function(inputResults,
                         InteractionCoeffcutoff=0.5){

  if(class(inputResults) != "IntLimResults") {
    stop("input data is not a IntLim class")
  }

  gene_metabolite_format_results = melt(inputResults@interaction.coefficients)
  colnames(gene_metabolite_format_results) = c("gene", "metabolite", "interaction")
  tofilter_sorted <- gene_metabolite_format_results[order(gene_metabolite_format_results$interaction),]

  filtered = tofilter_sorted[tofilter_sorted$interaction>InteractionCoeffcutoff | tofilter_sorted$interaction < -InteractionCoeffcutoff,]

  inputResults@filt.results = filtered
  return(inputResults)

}

#' Create results table, which includes significant gene:metabolite pairs, associated p-values,
#' and correlations in each category evaluated.
#'
#' @param inputResults IntLimResults object with model results (output of ProcessResults())
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata)
#' mytable <- CreateResultsTable(myres)
#' }
#' @export
#   CreateResultsTable <- function(inputResults) {
#        a<-inputResults@corr
#        a$cordiff<-round(abs(a[,3]-a[,4]),3)
#        a[,3]<-round(a[,3],2)
#        a[,4]<-round(a[,4],2)
#        p <- padj <- c()
#        if(inputResults@outcome=="metabolite") {
#                for (i in 1:nrow(a)) {
#                        g <- which(rownames(inputResults@interaction.pvalues) == a$gene[i])
#                        m <- which(colnames(inputResults@interaction.pvalues) == a$metab[i])
#                        if(length(g)==0 || length(m)==0) {p<-c(p,NA);padj<-c(padj,NA)} else {
#                                p <- c(p,inputResults@interaction.pvalues[g,m])
#				padj <- c(padj,inputResults@interaction.adj.pvalues[g,m])
#                              padj <- c(padj,inputResults@interaction.adj.pvalues[a$gene[i],a$metab[i]])
#                        }
#                }
#        } else if (inputResults@outcome=="gene") {
#               for (i in 1:nrow(a)) {
#                        g <- which(rownames(inputResults@interaction.pvalues) == a$gene[i])
#                        m <- which(colnames(inputResults@interaction.pvalues) == a$metab[i])
#                        if(length(g)==0 || length(m)==0) {p<-c(p,NA)} else {
 #                               p <- c(p,inputResults@interaction.pvalues[a$metab[i],a$gene[i]])
#				p <- c(p,inputResults@interaction.pvalues[g,m])
#                                padj <- c(padj,inputResults@interaction.adj.pvalues[m,g])
#                        }
#                }
#        }
#        else {stop("Outcome should be either 'metabolite' or 'gene'")}
#        a$pval <- p
#        a$adjpval <- padj
#        table<-a[order(a$adjpval,decreasing = TRUE),]
#	rownames(table) <- NULL
#        return(table)
#   }
#
#
