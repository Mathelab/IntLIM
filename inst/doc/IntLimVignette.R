## ----eval=FALSE----------------------------------------------------------
#  	runIntLimApp()

## ------------------------------------------------------------------------
     dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
     csvfile <- file.path(dir, "NCItestinput.csv")
     csvfile

## ------------------------------------------------------------------------
inputData <- IntLIM::ReadData(csvfile,metabid='id',geneid='id')
IntLIM::ShowStats(inputData)

## ------------------------------------------------------------------------
inputDatafilt <- IntLIM::FilterData(inputData,geneperc = 0.15)
IntLIM::ShowStats(inputDatafilt)

## ------------------------------------------------------------------------
IntLIM::PlotDistributions(inputData)

## ------------------------------------------------------------------------
IntLIM::PlotPCA(inputData,stype = "PBO_vs_Leukemia")

## ------------------------------------------------------------------------
myres <- IntLIM::RunIntLim(inputData,stype="PBO_vs_Leukemia")
IntLIM::DistPvalues(myres)

## ------------------------------------------------------------------------
myres <- IntLIM::ProcessResults(myres,inputData)
IntLIM::CorrHeatmap(myres)

## ------------------------------------------------------------------------
IntLIM::PlotGMPair(inputData,stype="PBO_vs_Leukemia","DLG4","(p-Hydroxyphenyl)lactic acid")

## ------------------------------------------------------------------------
sessionInfo()

