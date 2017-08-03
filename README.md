# IntLIM:  Integration through LInear Modeling
[![Build Status](https://travis-ci.org/Mathelab/IntLIM.svg?branch=master)](https://travis-ci.org/Mathelab/IntLIM)

## IntLIM

The goal of the IntLIM R package is to identify gene:metabolite relationships that are specific to a given phenotype (e.g. cancer vs non-cancer). For example, a given gene:metabolite pair could show a strong correlation in one phenotype (e.g. cancer) and no correlation in the other (e.g. non-cancer).  Users are expected to provide normalized gene expression and metabolite abundance data, as well as associated meta-information on the samples (at the minimum, users need to provide a phenotype column).  Currently, IntLIM requires the phenotype of interest to have 2 categories.  Optionally, users can also input meta-information on metabolites and genes (e.g. names, pathways).  

An example data set is provided within the package, and is a subset of the NCI-60 gene expression and metabolomics data (https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data).  

## IntLIM prerequisites
IntLIM is an R package and can be run on version >= 3.2.0. 

Download (or upgrade) R here: https://cloud.r-project.org/

RStudio (an interface to R than can make R easier to use) can be download here (not required): https://www.rstudio.com/products/rstudio/download3/

## Contact

If you encounter any problems running on the software, or find installation problems or bugs, please start an issue on the Issues tab or email Ewy Mathe at Ewy.Mathe@osumc.edu or Jalal Siddiqui at jalal.siddiqui@osumc.edu.  We are also very open to any comments, including how we can ameliorate the package.

## Installation from Github

To install IntLIM, simply type the following in the R terminal:
```
install.packages("devtools")
devtools::install_github("mathelab/IntLIM")
```

## Running IntLIM's user-friendly web app:

The package functions can be run direclty in the R console.  
Alternatively, to launch the web app, type th following in your R console:

```
library(IntLIM)
runIntLIMApp()
```
## Vignette
A detailed vignette can be found here:
https://mathelab.github.io/IntLIM/vignette.html

## Highcharts

All plots from this package use Highcharts: Highcharts (www.highcharts.com) is a Highsoft software product which is not free for commercial and Governmental use.
