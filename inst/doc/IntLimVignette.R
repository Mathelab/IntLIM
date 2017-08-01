## ----eval=FALSE----------------------------------------------------------
#  	runIntLimApp()

## ------------------------------------------------------------------------
     dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
     csvfile <- file.path(dir, "NCItestinput.csv")
     csvfile

## ------------------------------------------------------------------------
inputData <- IntLim::ReadData(csvfile,metabid='id',geneid='id')
IntLim::ShowStats(inputData)

## ------------------------------------------------------------------------
inputDatafilt <- IntLim::FilterData(inputData,geneperc = 0.15)
IntLim::ShowStats(inputDatafilt)

## ------------------------------------------------------------------------
IntLim::PlotDistributions(inputData)

## ------------------------------------------------------------------------
IntLim::PlotPCA(inputData,stype = "PBO_vs_Leukemia")

## ------------------------------------------------------------------------
myres <- IntLim::RunIntLim(inputData,stype="PBO_vs_Leukemia")
IntLim::DistPvalues(myres)

## ------------------------------------------------------------------------
myres2 <- IntLim::ProcessResults(myres,inputData)
IntLim::CorrHeatmap(myres2)

## ------------------------------------------------------------------------
IntLim::PlotGMPair(inputData,stype="PBO_vs_Leukemia","DLG4","(p-Hydroxyphenyl)lactic acid")

## ------------------------------------------------------------------------
sessionInfo()

