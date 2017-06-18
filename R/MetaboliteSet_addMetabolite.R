methods::setGeneric(
	"add_metabolite", 
	function(object, metabSet, warnings = TRUE, ...)
    	base::standardGeneric("add_metabolite")
)


methods::setMethod(
       f = "add_metabolite",
       signature = c("MultiDataSet", "MetaboliteSet"),
       definition = function(object, metabSet, warnings = TRUE, ...) {
           ## Add given MetaboliteSet as 'metabolite'
           object <- MultiDataSet::add_eset(object, metabSet, dataset.type = "metabolite", GRanges = NA, ...)
           return(object)
})
