methods::setMethod(
	f ="hist", 
	signature = "IntLimResults",
        definition = function(x, xlab='Unadjusted Interaction P-Values', main=NULL,...) {
              hist(x@interaction.pvalues, xlab=xlab, breaks=length(x@interaction.pvalues)/15, 
		main=main, ...)
          })

