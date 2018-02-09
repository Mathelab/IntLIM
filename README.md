# IntLIM:  Integration through LInear Modeling
[![Build Status](https://travis-ci.org/Mathelab/IntLIM.svg?branch=master)](https://travis-ci.org/Mathelab/IntLIM)
[![Build status](https://ci.appveyor.com/api/projects/status/1y05oo8y4v7r28bf?svg=true)](https://ci.appveyor.com/project/Mathelab/IntLIM/branch/master)

## IntLIM

Interpretation of metabolomics data is very challenging.  Yet it can be eased through integration of metabolmoics with other ‘omics’ data. The IntLIM package, which includes a user-friendly RShiny web app, aims to integrate metabolomics data with transcriptomic data.  Unlike other approaches, IntLIM is focused on understanding how specific gene-metabolite associations are affected by phenotypic features.  To this end, we develop a linear modeling approach that describes how gene-metabolite associations are affected by phenotype.  The workflow involves the following steps: 1) input gene expression/metabolomics data files, 2) filter data sets by gene and metabolite abundances and imputed values, 3) run the linear model to extract FDR-adjusted interaction p-values, 4) filter results by p-values and Spearman correlation differences, and 5) plot/visualize specific gene-metabolite associations. 

An example data set is provided within the package, and is a subset of the NCI-60 gene expression and metabolomics data (https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data).  The vignette outlines how to run the workflow.  

## IntLIM prerequisites

IntLIM is an R package and can be run on version >= 3.2.0. 

Download (or upgrade) R here: https://cloud.r-project.org/

RStudio (an interface to R than can make R easier to use) can be download here (not required): https://www.rstudio.com/products/rstudio/download3/

## Installation from Github

Prior to installing IntLIM, it is necessary to have the Bioconductor package *MultiDataSet* (Hernandez-Ferrer et al, 2017).  


The following command then installs *MultiDataSet*.

```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("MultiDataSet")
```

To install IntLIM, simply type the following in the R terminal:

```
install.packages("devtools")
library(devtools)
devtools::install_github("mathelab/IntLIM")
```
## Vignette

A detailed vignette can be found here:
https://mathelab.github.io/IntLIM/IntLIMVignette.html

## Formatted Data and Analysis Codes

Formatted data and codes to reproduce the NCI-60 analyses can be obtained from the following GitHub repository:

https://github.com/Mathelab/NCI60_GeneMetabolite_Data

Formatted data and codes to reproduce the breast cancer analyses can be obtained from the following GitHub repository:

https://github.com/Mathelab/BreastCancerAmbs_GeneMetabolite_Data

## Running IntLIM's user-friendly web app:

The package functions can be run directly in the R console.  
Alternatively, to launch the web app, type the following in your R console:

```
library(IntLIM)
runIntLIMApp()
```

## Contact

If you encounter any problems running on the software, or find installation problems or bugs, please start an issue on the Issues tab or email Ewy Mathe at Ewy.Mathe@osumc.edu or Jalal Siddiqui at jalal.siddiqui@osumc.edu.  We are also very open to any comments, including how we can improve and ameliorate the package.
