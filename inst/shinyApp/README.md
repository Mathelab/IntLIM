# Welcome to the IntLIM shiny app!

The goal of this app is to provide users with a user-friendly platform for integrating gene expression and metabolomics data.  Specifically, the software finds gene:metabolite relationships that are specific to a given phenotype (e.g. cancer vs non-cancer). For example, a given gene:metabolite pair could show a strong correlation in one phenotype (e.g. cancer) and no correlation in the other (e.g. non-cancer). 

More details can be found in <a href="https://arxiv.org/abs/1802.10588" target="_blank">our manuscript</a>.

## Getting started (loading in data)

__*Please be sure that all files noted in the CSV file, including the CSV file, are in the same folder. Do not include path names in the filenames.*__

Users need to input a CSV file with two required columns: 'type' and 'filenames'.
The CSV file is expecected to have the following 2 columns and 6 rows:
1. type,filenames
2. metabData,myfilename
3. geneData,myfilename
4. metabMetaData,myfilename (optional)
5. geneMetaData,myfilename (optional)
6. sampleMetaData,myfilename"

Note also that the input data files should be in a specific format:
- metabData: rows are metabolites, columns are samples; the first row is assumed to have sample ids and these ids should be unique; the first column is assumed to have feature ids and those should be unique.
- geneData: rows are genes, columns are samples; the first row is assumed to have sample ids and these ids should be unique; the first column is assumed to have feature ids and those should be unique.
- metabMetaData: rows are metabolites, features are columns
- geneMetaData: rows are genes, features are columns
- sampleMetaData: rows are samples, features are columns

*NOTE*: The first column of the sampleMetaData file is assumed to be the sample id, and those sample ids should match the *first row* of metabData and geneData (e.g. it is required that all sample ids in the metabData and geneData are also in the sampleMetaDatafile).

## Test data
The package includes a reduced set of the original NCI-60 dataset.  The CSV input file location for this test dataset can be located by typing the following in the R console:
```
     dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
     csvfile <- file.path(dir, "NCItestinput.csv")
     csvfile
```
Please see the vignette at [https://mathelab.github.io/IntLIM/IntLIMVignette.html](https://mathelab.github.io/IntLIM/IntLIMVignette.html) for additional information.

In addition, the original datasets for the NCI-60 and the breast cancer dataset used in the manuscript can be found HEREXXX.

## Contact

If you have any questions, comments, or concerns on how to use IntLIM please contact Ewy Mathe at ewy.mathe@osumc.edu or  Jalal Siddiqui at jalal.siddiqui@osumc.edu.
