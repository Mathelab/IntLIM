methods::setGeneric(
	"getStats", 
	function(IntLimObject, ...)
    	base::standardGeneric("getStats")
)


methods::setMethod(
       f = "getStats",
       signature = c("MultiDataSet"),
       definition = function(IntLimObject, ...) {
		incommon<-MultiDataSet::commonSamples(IntLimObject)
	       mystats <- data.frame(Num_Genes = nrow(Biobase::fData(IntLimObject[["expression"]])),
		Num_Metabolites = nrow(Biobase::fData(IntLimObject[["metabolite"]])),
		Num_Samples_withExpression = ncol(Biobase::assayDataElement(IntLimObject[["expression"]], 
			'exprs')),
		Num_Samples_withExpression = ncol(Biobase::assayDataElement(IntLimObject[["metabolite"]], 
			'metabData')),
		Num_Samples_inCommon = ncol(Biobase::assayDataElement(incommon[["expression"]], 'exprs'))
		)
           return(mystats)
})
