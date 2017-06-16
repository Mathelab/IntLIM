#' MetaboliteSet class extenstion for MultiDataSet object
#'
#' This class is specific to metabolomics data and will then be stored into the class
#' \code{MultiDataSet}, which stores multiple dataset types (e.g. metabolite and gene levels
#' in this case.
#'
#' @name MetaboliteSet-class
#' @rdname MetaboliteSet-class
#' @exportClass MetaboliteSet
#' @slot eSet List of eSet elements

methods::setClass (
	Class = "MetaboliteSet",
	contains = "eSet",
	prototype = methods::prototype(methods::new("VersionedBiobase",
		versions = c(Biobase::classVersion("eSet"),
		MetaboliteSet = "1.0.0")))
           # where=topenv(parent.frame())
)
