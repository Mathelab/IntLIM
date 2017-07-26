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
				pvalcutoff=0.1,
				diffcorr=0.5,
				corrtype="spearman"){

	if(inputResults@outcome == "metabolite") {
		mydat <- reshape2::melt(inputResults@interaction.adj.pvalues)}
	else if (inputResults@outcome == "gene") {
                mydat <- reshape2::melt(t(inputResults@interaction.adj.pvalues))}

	keepers <- which(mydat$value <= pvalcutoff)

	incommon <- getCommon(inputData,inputResults@stype)
	p <- incommon$p
	gene <- incommon$gene
	metab <- incommon$metab
	if(length(unique(p)) !=2) {
 		stop(paste("IntLim currently requires only two categories.  Make sure the column",inputResults@stype,"only has two unique values"))
    }

	gp1 <- which(p == unique(p)[1])
	cor1 <- as.numeric(apply(mydat[keepers,],1,function(x) {
		stats::cor(as.numeric(gene[as.character(unlist(x[1])),gp1]),
			as.numeric(metab[as.character(unlist(x[2])),gp1]),method=corrtype)}))

	gp2 <- which(p == unique(p)[2])
        cor2 <- as.numeric(apply(mydat[keepers,],1,function(x) {
	         stats::cor(as.numeric(gene[as.character(unlist(x[1])),gp2]),
                        as.numeric(metab[as.character(unlist(x[2])),gp2]),method=corrtype)}))

        mydiffcor = abs(cor1-cor2)

	keepers2 <- which(mydiffcor > diffcorr)

	inputResults@corr <- data.frame(metab=as.character(mydat[keepers[keepers2],2]), 
		gene=as.character(mydat[keepers[keepers2],1]))
	inputResults@corr <- cbind(inputResults@corr,cor1[keepers2],cor2[keepers2])
	colnames(inputResults@corr)[3:4]=setdiff(as.character(unlist(unique(p))),"")

return(inputResults)
}
