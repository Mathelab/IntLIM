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

if(!require("Biobase")) install.packages("Biobase") 

methods::setClass (
	Class = "MetaboliteSet",
	contains = "eSet",
	#prototype = methods::prototype(methods::new("VersionedBiobase",
#		versions = c(Biobase::classVersion("eSet"),
#		MetaboliteSet = "1.0.0")))
           # where=topenv(parent.frame())
)

#' MetaboliteSet class extenstion for MultiDataSet object
#'
#' This class is specific to metabolomics data and will then be stored into the class
#' \code{MultiDataSet}, which stores multiple dataset types (e.g. metabolite and gene levels
#' in this case.
#'
#' @name IntLimModel-class
#' @rdname IntLimModel-class
#' @exportClass IntLimModel
#' @slot call function call
#' @slot model model formula
#' @slot F.statistics numeric
#' @slot p.values numeric
#' @slot coefficients matrix
#' @slot predictions matrix
#' @slot sse numeric
#' @slot ssr numeric
#' @slot df numeric
#' @slot std.error.coeff matrix
#' @slot t.value.coeff matrix
#' @slot p.value.coeff
methods::setClass('IntLimModel',
         representation(call='call',
                        model='formula',
                        F.statistics='numeric',
                        p.values='numeric',
                        coefficients='matrix',
                        predictions='matrix',
                        sse='numeric',
                        ssr='numeric',
                        df='numeric', 
                        std.error.coeff = 'matrix',
                        t.value.coeff = 'matrix',
                        p.value.coeff = 'matrix'))


