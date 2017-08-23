#' Get some stats after reading in data
#'
#' @include MultiDataSet_extendedfunctions.R
#'
#' @param IntLimObject output of ReadData()
#' @return data.frame with some # of samples, features, etc.
#'
#' @examples
#' dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' ShowStats(mydata)
#' @export
ShowStats <- function(IntLimObject) {
        return(getStats(IntLimObject))
}
