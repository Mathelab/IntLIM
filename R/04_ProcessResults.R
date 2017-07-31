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
#' @return IntResults object with model results (now includes correlations)
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata)
#' }
#' @export
ProcessResults <- function(inputResults,
				inputData,
				pvalcutoff=0.05,
				diffcorr=0.5,
				corrtype="spearman"){

	if(inputResults@outcome == "metabolite") {
		mydat <- reshape2::melt(inputResults@interaction.adj.pvalues)}
	else if (inputResults@outcome == "gene") {
                mydat <- reshape2::melt(t(inputResults@interaction.adj.pvalues))}

	keepers <- which(mydat$value <= pvalcutoff)
	print(length(keepers))

	incommon <- getCommon(inputData,inputResults@stype)
	p <- incommon$p
	gene <- incommon$gene
	metab <- incommon$metab
	if(length(unique(p)) !=2) {
 		stop(paste("IntLim currently requires only two categories.  Make sure the column",inputResults@stype,"only has two unique values"))
    }

	print("Processing gp1")
	if(pvalcutoff == 1)
	gp1 <- which(p == unique(p)[1])
	cor1.m <- cor(t(gene),t(metab),method=corrtype)
	if(pvalcutoff == 1) {temp <- reshape::melt(cor1.m); fincor1 <- temp$value
		} else {
		fincor1 <- as.numeric(lapply(keepers,function(x) 
			cor1.m[as.character(mydat$Var1[x]),as.character(mydat$Var2[x])]))
		}

#	cor1 <- reshape::melt(cor1.m)
#	myind <- as.numeric(lapply(keepers,function(x) 
#		intersect(which(as.character(cor1$X1)==as.character(mydat$Var1[x])),
#			which(as.character(cor1$X2)==as.character(mydat$Var2[x])))))
#	fincor1 <- cor1[myind,]

	print("Processing gp2")
	gp2 <- which(p == unique(p)[1])
        cor2.m <- cor(t(gene),t(metab),method=corrtype)
        if(pvalcutoff == 1) {temp <- reshape::melt(cor2.m); fincor2 <- temp$value
                } else {
			fincor2 <- as.numeric(lapply(keepers,function(x)
                	cor2.m[as.character(mydat$Var1[x]),as.character(mydat$Var2[x])]))
	}

	#gp1 <- which(p == unique(p)[1])
	#cor1 <- as.numeric(apply(mydat[keepers,],1,function(x) {
	#	stats::cor(as.numeric(gene[as.character(unlist(x[1])),gp1]),
	#		as.numeric(metab[as.character(unlist(x[2])),gp1]),method=corrtype)}))

	#gp2 <- which(p == unique(p)[2])
        #cor2 <- as.numeric(apply(mydat[keepers,],1,function(x) {
	#         stats::cor(as.numeric(gene[as.character(unlist(x[1])),gp2]),
        #                as.numeric(metab[as.character(unlist(x[2])),gp2]),method=corrtype)}))


        mydiffcor = abs(fincor1-fincor2)

	keepers2 <- which(mydiffcor > diffcorr)

	inputResults@corr <- data.frame(metab=as.character(mydat[keepers[keepers2],2]), 
		gene=as.character(mydat[keepers[keepers2],1]))
	inputResults@corr <- cbind(inputResults@corr,fincor1[keepers2],fincor2[keepers2])
	colnames(inputResults@corr)[3:4]=setdiff(as.character(unlist(unique(p))),"")

return(inputResults)
}


#' Create results table, which includes significant gene:metabolite pairs, associated p-values, 
#' and correlations in each category evaluated.
#'
#' @param inputResults IntLimResults object with model results (output of ProcessResults())
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata)
#' mytable <- CreateResultsTable(myres)
#' }
#' @export
   CreateResultsTable <- function(inputResults) {
        a<-inputResults@corr
        a$cordiff<-round(abs(a[,3]-a[,4]),3)
        a[,3]<-round(a[,3],2)
        a[,4]<-round(a[,4],2)
        p <- padj <- c()
        if(inputResults@outcome=="metabolite") {
                for (i in 1:nrow(a)) {
                        g <- which(rownames(inputResults@interaction.pvalues) == a$gene[i])
                        m <- which(colnames(inputResults@interaction.pvalues) == a$metab[i])
                        if(length(g)==0 || length(m)==0) {p<-c(p,NA);padj<-c(padj,NA)} else {
                                p <- c(p,inputResults@interaction.pvalues[g,m])
				padj <- c(padj,inputResults@interaction.adj.pvalues[g,m]) 
#                              padj <- c(padj,inputResults@interaction.adj.pvalues[a$gene[i],a$metab[i]])
                        }
                }
        } else if (inputResults@outcome=="gene") {
               for (i in 1:nrow(a)) {
                        g <- which(rownames(inputResults@interaction.pvalues) == a$gene[i])
                        m <- which(colnames(inputResults@interaction.pvalues) == a$metab[i])
                        if(length(g)==0 || length(m)==0) {p<-c(p,NA)} else {
 #                               p <- c(p,inputResults@interaction.pvalues[a$metab[i],a$gene[i]])
				p <- c(p,inputResults@interaction.pvalues[g,m])
                                padj <- c(padj,inputResults@interaction.adj.pvalues[m,g])
                        }
                }
        }
        else {stop("Outcome should be either 'metabolite' or 'gene'")}
        a$pval <- p
        a$adjpval <- padj
        table<-a[order(a$adjpval,decreasing = TRUE),]
	rownames(table) <- NULL
        return(table)
   }


