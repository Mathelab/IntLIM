## ----eval=FALSE----------------------------------------------------------
#  	runIntLimApp()

## ------------------------------------------------------------------------
     dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
     csvfile <- file.path(dir, "NCItestinput.csv")
     csvfile

## ------------------------------------------------------------------------
inputData <- IntLim::ReadData(csvfile,metabid='id',geneid='id')
IntLim::OutputStats(inputData)

## ------------------------------------------------------------------------
inputDatafilt <- IntLim::FilterData(inputData,geneperc=15)
IntLim::OutputStats(inputDatafilt)

## ------------------------------------------------------------------------
IntLim::PlotDistributions(inputData)

## ------------------------------------------------------------------------
myres <- IntLim::RunIntLim(inputData,stype="PBO_vs_Leukemia")
IntLim::DistPvalues(myres)

## ------------------------------------------------------------------------
myres <- IntLim::ProcessResults(myres,inputData)
IntLim::CorrHeatmap(myres)

## ------------------------------------------------------------------------
sessionInfo()

