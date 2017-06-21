#' Get some stats after reading in data
#'
#' @include MultiDataSet_extendedfunctions.R
#'
#' @param IntLimObject output of ReadData()
#' @return data.frame with some # of samples, features, etc.
#'
#' @examples
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "test.csv")
#' mydata <- ReadData(csvfile,metabid='BIOCHEMICAL',geneid='X')
#' OutputStats(mydata)
#' @export
OutputStats <- function(IntLimObject) {
        return(getStats(IntLimObject))
}
