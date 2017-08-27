## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval = FALSE--------------------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("MultiDataSet")

## ----eval = FALSE--------------------------------------------------------
#  install.packages(devtools)
#  install_github(“/mathelab/IntLIM”)

## ------------------------------------------------------------------------
library(IntLIM)

## ------------------------------------------------------------------------
     dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
     csvfile <- file.path(dir, "NCItestinput.csv")
     csvfile

## ------------------------------------------------------------------------
inputData <- IntLIM::ReadData(inputFile = csvfile,metabid='id',geneid='id')

## ------------------------------------------------------------------------
IntLIM::ShowStats(IntLimObject = inputData)

## ------------------------------------------------------------------------
inputDatafilt <- IntLIM::FilterData(inputData,geneperc = 0.10, metabmiss = 0.80)
IntLIM::ShowStats(inputDatafilt)

## ------------------------------------------------------------------------
IntLIM::PlotDistributions(inputData = inputDatafilt)

## ------------------------------------------------------------------------
IntLIM::PlotPCA(inputData = inputDatafilt,stype = "PBO_vs_Leukemia")

## ------------------------------------------------------------------------
myres <- IntLIM::RunIntLim(inputData = inputDatafilt,stype="PBO_vs_Leukemia")

## ------------------------------------------------------------------------
IntLIM::DistPvalues(IntLimResults = myres)

## ------------------------------------------------------------------------
IntLIM::pvalCorrVolcano(inputResults = myres, inputData = inputDatafilt, diffcorr = 0.5, pvalcutoff = 0.05)

## ------------------------------------------------------------------------
myres <- IntLIM::ProcessResults(inputResults = myres, inputData = inputDatafilt, pvalcutoff = 0.10, diffcorr = 0.5)

## ------------------------------------------------------------------------
IntLIM::CorrHeatmap(myres)

## ------------------------------------------------------------------------
IntLIM::PlotGMPair(inputDatafilt,stype="PBO_vs_Leukemia","DLG4","(p-Hydroxyphenyl)lactic acid")

## ----eval = FALSE--------------------------------------------------------
#  IntLIM::OutputData(inputData=inputDatafilt,filename="~/FilteredData.zip")
#  OutputResults(inputResults=myres,filename="~/MyResults.zip")
#  

## ----eval=FALSE----------------------------------------------------------
#  	runIntLIMApp()

