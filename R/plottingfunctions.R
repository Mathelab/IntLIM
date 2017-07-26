#' Get some stats after reading in data
#'
#' @import magrittr
#'
#' @include MultiDataSet_extendedfunctions.R
#'
#' @param inputData IntLimObject output of ReadData()
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @return a highcharter object
#'
#' @examples
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' PlotDistributions(mydata)
#' @export
PlotDistributions <- function(inputData,viewer=T,
#	palette = c("#C71585", "#00E5EE")) {
        palette="Set1"){

      if (length(palette) == 2) {
        cols <- c(palette)
      }
      else if (length(palette) == 1) {
        cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
      }
      else {
        stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
      }

    categ <- c("Genes","Metabolites")

	mygene <- as.data.frame(Biobase::assayDataElement(inputData[["expression"]],'exprs'))
	toplot <- suppressMessages(reshape2::melt(mygene))
        df <- dplyr::data_frame(value = toplot$value, by = toplot$variable) %>% dplyr::group_by_("by") %>%
 	       dplyr::do(data = grDevices::boxplot.stats(.$value))
#	names(df$data) <- df$by
#        df$color <- df$by
	bxps <- purrr::map(df$data, "stats")
	outs <- purrr::map2_df(seq(nrow(df)), df$data, function(x, y) {
            if (length(y$out) > 0)
                d <- dplyr::data_frame(x = x - 1, y = y$out)
            else d <- dplyr::data_frame()
            d
        })
# To try to get the gene names of outliers, would have to go back and get the gene names from original data frame and put htem in outs$color
	boxplotOptions <- list(
          fillColor = '#ffffff',
          lineWidth = 2,
          medianColor = '#000000',
          medianWidth = 2,
          stemColor = '#000000',
          stemDashStyle = 'dot',
          stemWidth = 1,
          whiskerColor = '#000000',
          whiskerLength = '20%',
          whiskerWidth = 3)

	g <- highcharter::highchart(width = 750, height = 750 ) %>%
      highcharter::hc_title(text = "Gene Expression",
               style = list(color = '#2E1717',
			fontWeight = 'bold', fontSize = "20px")) %>%
      highcharter::hc_plotOptions(
        boxplot = boxplotOptions
        ) %>%
      hc_add_series(data = bxps,type="boxplot",color=cols[1],showInLegend=FALSE) %>%
      highcharter::hc_add_series(data=list_parse(outs),type="scatter",color=cols[1],showInLegend=FALSE) %>%
#		name = str_trim(paste(list(...)$name, "outliers")),
#                type = "scatter") #marker = list(...)) %>%
#		tooltip = list(headerFormat = "<span>{point.key}</span><br/>")) %>%
#                  pointFormat = "<span style='color:{point.color}'></span> \nOutlier: <b>{point.y}</b><br/>")) %>%
      highcharter::hc_yAxis(title = list(text = "log(expression)",
                            style = list(fontSize = "13px")),
               labels = list(format = "{value}")) %>%
      highcharter::hc_xAxis(labels="") %>%
      highcharter::hc_tooltip(valueDecimals = 2) %>%
      highcharter::hc_exporting(enabled = TRUE)

	mymetab <- Biobase::assayDataElement(inputData[["metabolite"]],'metabData')
	toplot <- suppressMessages(reshape2::melt(mymetab))
        df <- dplyr::data_frame(value = toplot$value, by = toplot$variable) %>% 
		dplyr::group_by_("by") %>%
               dplyr::do(data = grDevices::boxplot.stats(.$value))
        bxps <- purrr::map(df$data, "stats")
        outs <- purrr::map2_df(seq(nrow(df)), df$data, function(x, y) {
            if (length(y$out) > 0)
                d <- dplyr::data_frame(x = x - 1, y = y$out)
            else d <- dplyr::data_frame()
            d
        })

        m <- highcharter::highchart(width = 750, height = 750 ) %>%
      highcharter::hc_title(text = "Metabolite Levels",
               style = list(color = '#2E1717',
                            fontWeight = 'bold', fontSize = "20px")) %>%
      highcharter::hc_plotOptions(
        boxplot = boxplotOptions
        ) %>%
      highcharter::hc_add_series(data = bxps,type="boxplot",color=cols[2],showInLegend=FALSE) %>%
      highcharter::hc_add_series(data=list_parse(outs),type="scatter",color=cols[2],showInLegend=FALSE) %>%
      
      highcharter::hc_yAxis(title = list(text = "log(abundances)",
                            style = list(fontSize = "13px")),
               labels = list(format = "{value}")) %>%
      highcharter::hc_xAxis(labels="") %>%
      highcharter::hc_tooltip(valueDecimals = 2) %>%
      highcharter::hc_exporting(enabled = TRUE)

  if (viewer == TRUE) {
    p <-
      htmltools::browsable(highcharter::hw_grid(g, m, ncol = 2, rowheight = 550))
  }
  else {
    p <- highcharter::hw_grid(g, m)
  }
  return(p)
}

#' PCA plots of data for QC
#'
#' @import magrittr
#'
#' @include MultiDataSet_extendedfunctions.R
#'
#' @param inputData IntLimObject output of ReadData()
#' @param stype category to color-code by (can be more than two categories)
#' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
#' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param common whether or not samples that are in common between the metabolite and gene expression datasets should be plotted (T/F); default is TRUE
#' @return a highcharter object
#'
#' @examples
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' PlotPCA(mydata,stype = "PBO_vs_Leukemia")
#' @export
PlotPCA <- function(inputData,viewer=T,stype=NULL,common=T,
        palette = "Set1") {

    categ <- c("Genes","Metabolites")
    
        if(is.null(stype)) {
		warning("The resulting PCA plot is not color-coded because you did not provide a category in 'stype'")
		mytype <- NULL
        } else if (length(intersect(colnames(Biobase::pData(inputData[["metabolite"]])),stype))!=1) {
		stop(paste0("You provided ",stype, "as your stype variable but it does not exist in your data"))
        } else {
        	mytype <- as.character(Biobase::pData(inputData[["metabolite"]])[,stype])
                numcateg <- length(unique(mytype))
                if(length(palette) >= 2) {
                           cols <- palette 
                } else {
                if(numcateg == 1) {
                       if(length(palette)==1) {cols <- RColorBrewer::brewer.pal(3, palette)[1]
                       } else {stop("palette should be an RColorBrewer palette or a vector of colors")}
                } else if (numcateg == 2) {
                      if(length(palette)==1) {cols <- RColorBrewer::brewer.pal(numcateg, palette)[1:2]
                      } else {stop("palette should be an RColorBrewer palette or a vector of colors")}
                } else if (numcateg > 2) {
                      if(length(palette)==1) {cols <- RColorBrewer::brewer.pal(numcateg, palette)
                      } else {stop("palette should be an RColorBrewer palette or a vector of colors")}
                } else {stop("There are no values in your 'stype' column")}
               }
        }      


	if(common==T) {
		if(is.null(stype)) {
			incommon <- getCommon(inputData)
			mygene <- incommon$gene
			gpca <- stats::prcomp(t(mygene),center=F,scale=T)
			mymetab <- incommon$metab
			mpca <- stats::prcomp(t(mymetab),center=F,scale=T)
			gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),color=rep("blue",nrow(gpca$x)))
			mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),color=rep("blue",nrow(mpca$x)))
		} else {
			incommon <- getCommon(inputData,stype)
			mygene <- incommon$gene
			mymetab <- incommon$metab
			mytype <- incommon$p
			uniqtypes <- unique(mytype)
			mycols <- as.character(mytype)	
			for (i in 1:numcateg) {
				mycols[which(mytype==uniqtypes[i])] <- cols[i]
			}
			gpca <- stats::prcomp(t(mygene),center=F,scale=T)
			mpca <- stats::prcomp(t(mymetab),center=F,scale=T)
			gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),label=mytype,color=mycols)
			mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),label=mytype,color=mycols)
		}
	} else { # common == F
		mygene <- as.data.frame(Biobase::assayDataElement(inputData[["expression"]],'exprs'))
		mymetab <- Biobase::assayDataElement(inputData[["metabolite"]],'metabData')
        	gpca <- stats::prcomp(t(mygene),center=F,scale=T)
		mpca <- stats::prcomp(t(mymetab),center=F,scale=T)
#		percvar=round((gpca$sdev)^2 / sum(gpca$sdev^2)*100,2)
		if(!is.null(stype)) {
			gtypes <- as.character(Biobase::pData(inputData[["expression"]])[,stype])
			mtypes <- as.character(Biobase::pData(inputData[["metabolite"]])[,stype])
        		uniqtypes <- unique(c(mtypes,gtypes))
        		gcols <- as.character(gtypes)
			mcols <- as.character(mtypes)
        		for (i in 1:numcateg) {
				gcols[which(gtypes==uniqtypes[i])] <- cols[i]
				mcols[which(mtypes==uniqtypes[i])] <- cols[i]
        		}
		        # Deal with missing values or ""
		        if(length(which(gtypes==""))>0) {
				gcols[which(gtypes=="")]="grey"
        		        gtypes[which(gtypes=="")]="NA"
        		}
			if (length(which(mtypes==""))>0) {
                                mcols[which(mtypes=="")]="grey"
                                mtypes[which(mtypes=="")]="NA"
			}
		        if(length(which(is.na(gtypes)))>0) {
                		gcols[which(is.na(gtypes))]="grey"
			}        
                        if(length(which(is.na(mtypes)))>0) {
                                mcols[which(is.na(mtypes))]="grey"
                        }

			gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),label=gtypes,color=gcols)
                        mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),label=mtypes,color=mcols)
		} else { #stype is null
			gtoplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),color=rep("blue",nrow(gpca$x)))
                        mtoplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),color=rep("blue",nrow(mpca$x)))
		}
	} # end common == F

        mds <- list_parse(mtoplot)
	gds <- list_parse(gtoplot)
	mpercvar=round((mpca$sdev)^2 / sum(mpca$sdev^2)*100,2)
	gpercvar=round((gpca$sdev)^2 / sum(gpca$sdev^2)*100,2)

	pg <- highcharter::highchart(width = 350, height = 350 ) %>%
		highcharter::hc_title(text="PCA of genes") %>%
		highcharter::hc_xAxis(title=list(text=paste0("PC1:",round(gpercvar[1],1),"%"))) %>%
		highcharter::hc_yAxis(title=list(text=paste0("PC2:",round(gpercvar[2],2),"%"))) %>%
		hc_chart(zoomType = "xy") %>% 
		highcharter::hc_add_series(data=gds,type="scatter",col=cols[1],
			tooltip = list(headerFormat="", 
			  pointFormat=paste("{point.label}","{point.z}")),
			showInLegend=FALSE)
#		dataLabels= list(enabled = TRUE, format = "{point.label}"),

        pm <- highcharter::highchart(width = 350, height = 350 ) %>%
                highcharter::hc_title(text="PCA of metabolites") %>%
                highcharter::hc_xAxis(title=list(text=paste0("PC1:",round(mpercvar[1],1),"%"))) %>%
                highcharter::hc_yAxis(title=list(text=paste0("PC2:",round(mpercvar[2],2),"%"))) %>%
                hc_chart(zoomType = "xy") %>%
                highcharter::hc_add_series(data=mds,type="scatter",col=cols[1],
                        tooltip = list(headerFormat="",
                          pointFormat=paste("{point.label}","{point.z}")),
                        showInLegend=FALSE)


         if (viewer == TRUE) {
    		p <-
      		htmltools::browsable(highcharter::hw_grid(pg, pm, ncol = 2, rowheight = 550))
  	} else {
    		p <- highcharter::hw_grid(pg, pm)
  	}
  return(p)

}



#' Visualize the distribution of unadjusted p-values from linear models
#'
#' @include IntLimResults_extendedfunctions.R
#'
#' @param IntLimResults output of RunIntLim()
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' DistPvalues(myres)
#' }
#' @export
DistPvalues<- function(IntLimResults) {
    y<-as.numeric(IntLimResults@interaction.pvalues)
    hchart(y)
}


#' Plot correlation heatmap
#'
#' @import magrittr
#' @import highcharter
#'
#' @param inputResults IntLimResults object (output of ProcessResults())
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCItestinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' myres <- RunIntLim(mydata,stype="PBO_vs_Leukemia")
#' myres <- ProcessResults(myres,mydata)
#' CorrHeatmap(myres)
#' }
#' @export
CorrHeatmap <- function(inputResults,viewer=T) {
type <- cor <- c()

	if(nrow(inputResults@corr)==0) {
		stop("Make sure you run ProcessResults before making the heatmap")
	}
		temp <- inputResults@corr
		toplot <- data.frame(name=paste(temp[,1],temp[,2],sep=" vs "),
			temp[,3:4])
		suppressMessages(
			meltedtoplot <- tidyr::gather(
				toplot,
				type,cor,colnames(toplot)[2],colnames(toplot)[3]))

		#all possible values of X (type) and Y (name)
  		theXAxis <- as.character(meltedtoplot[, "type"])
		theYAxis <- as.character(meltedtoplot[, "name"])

		  #unique values of X and Y
		  theUniqueY <- as.character(unique(theYAxis))
		  theUniqueX <- as.character(unique(theXAxis))

		  # Substitute words with position on the meatrix
		  for (i in 1:length(theUniqueY)){
		    num <- which(theYAxis == theUniqueY[i])
		    theYAxis[num] <- i
		  }
		  for (i in 1:length(theUniqueX)) {
		    num <- which(theXAxis == theUniqueX[i])
		    theXAxis[num] <- i
		  }

		  #create final formatting
		  dataforHeatmap <- as.data.frame(cbind(
		    as.numeric(theXAxis),
		    as.numeric(theYAxis),
		    as.numeric(meltedtoplot$cor)
		#as.numeric(format(meltedtoplot$cor,scientific=T,digits=2))
		  ))

		  formattedHeatmapData <- list_parse2(dataforHeatmap)

		  fntltp <- JS(
		    "function(){
		    return 'cor='+this.point.value;
		    }")

		p <- highchart(width = 800, height = 700) %>%
		    hc_chart(type = "heatmap", spacingRight = 160) %>%
		    hc_title(text = "Correlation Heatmap",
		             style = list(color = '#2E1717',fontSize = '20px',
		                          fontWeight = 'bold')) %>%
		    hc_xAxis(categories = c("",unique(as.character(meltedtoplot[, "type"]))),
		             labels = list(style = list(fontSize = '10px'))) %>%
		    hc_yAxis(categories = c("",unique(as.character(meltedtoplot[, "name"]))), 
				labels = list(style = list(fontSize = '10px'))) %>%
		    hc_add_series(data = formattedHeatmapData) %>% 
		    hc_tooltip(formatter = fntltp, valueDecimals = 2) %>%
		    hc_colorAxis(stops = color_stops(2, colors = c("#5097D1", "#DEEFF5")),
		                 min = min(as.numeric(dataforHeatmap[ , 3]), na.rm = T),
		                 max = max(as.numeric(dataforHeatmap[ , 3]), na.rm = T)) %>%
		    hc_legend(
		      enabled = TRUE,
		      layout = "vertical",
		      align = "right",
		      verticalAlign = "top",
		      floating = FALSE,
		      maxWidth = 200,
		      x = -10, # 90
		      y = 100, # 70
		      padding = 2
		      #title = list(text="p-value")
		    ) %>%
	 	      hc_exporting(enabled = TRUE)
	return(p)	
}

#' scatter plot of gene-metabolite pairs (based on user selection)
#'
#' @import magrittr
#' @import highcharter
#'
#' @param inputData IntLimObject output of ReadData() or FilterData()
#' @param stype category to color-code by
##' @param palette choose an RColorBrewer palette ("Set1", "Set2", "Set3",
##' "Pastel1", "Pastel2", "Paired", etc.) or submit a vector of colors
#' @param geneName string of select geneName
#' @param viewer whether the plot should be displayed in the RStudio viewer (T) or
#' in Shiny/Knittr (F)
#' @param metabName string of select metabName
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCIinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' PlotGMPair(mydata,stype="PBO_vs_Leukemia","DLG4","(p-Hydroxyphenyl)lactic acid")
#' 
#' }
#' @export
PlotGMPair<- function(inputData,stype=NULL,geneName,metabName,palette = "Set1",
	viewer=T) {

      if(is.null(stype)) {
	stop("Users must define stype which defines the categories to be compared (e.g. tumor vs non-tumor).  This could be the same parameter that was used to run RunIntLim()")
	}
      if (length(palette) == 2) {
        cols <- c(palette)
      }
      else if (length(palette) == 1) {
        cols <- RColorBrewer::brewer.pal(3, palette)[1:2]
      }
      else {
        stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
      }
   
   if (class(inputData) != "MultiDataSet") {
        stop("input data is not a MultiDataSet class")
    }

    incommon <- getCommon(inputData,stype)

	if(is.null(stype)) {
                stop("A category to colorcode by (e.g. stype) must be provided")
        } else if (length(intersect(colnames(Biobase::pData(inputData[["metabolite"]])),stype))!=1) {
                stop(paste0("You provided ",stype, "as your stype variable but it does not exist in your data"))}
        else {
                mytypes <- incommon$p
        }

    gene<-incommon$gene
    if(length(which(rownames(gene)==geneName))>0) {
	    sGene<-gene[geneName,]
    } else {
	stop(paste0("The gene ",geneName," was not found in your data"))
    }
    
    metab<-incommon$metab
    if(length(which(rownames(metab)==metabName))>0) {
    	sMetab<-as.numeric(metab[metabName,])
    } else {
	stop(paste0("The metabolite ",metabName," was not found in your data"))
    }

    if(length(unique(mytypes))!=2) {
	stop(paste0("The group selected, '",stype,"', should only contain two different categories"))
    }   
 
    mycols <- as.character(mytypes)
    mycols[which(mytypes==unique(mytypes)[1])] <- cols[1]
    mycols[which(mytypes==unique(mytypes)[2])] <- cols[2]
    
    data<-data.frame(x=sGene,y=sMetab,z=colnames(gene),label=mytypes,color=mycols)

#    data<-data[data$label!="",]
    #data$type <- factor(data$type)

    max<- max(data$x)
    min<-min(data$x)

    m1<-stats::glm(data$y[which(data$label==mytypes[1])]~data$x[which(data$label==mytypes[1])])
    line1<-data.frame(x=c(max,min),
	y=c(as.numeric(m1$coefficients[2])*max+as.numeric(m1$coefficients[1]),
		as.numeric(m1$coefficients[2])*min+as.numeric(m1$coefficients[1])))
    m2<-stats::glm(data$y[which(data$label==mytypes[2])]~data$x[which(data$label==mytypes[2])])
    line2<-data.frame(x=c(max,min),
	y=c(as.numeric(m2$coefficients[2])*max+as.numeric(m2$coefficients[1]),
		as.numeric(m2$coefficients[2])*min+as.numeric(m2$coefficients[1])))

    ds <- list_parse(data)
    #cols=c("blue","pink")

        hc <- highcharter::highchart(width = 350, height = 350 ) %>%
                highcharter::hc_title(text="Gene:metabolite scatterplot") %>%
                highcharter::hc_xAxis(title=list(text=geneName)) %>%
                highcharter::hc_yAxis(title=list(text=metabName)) %>%
                hc_chart(zoomType = "xy") %>%
                highcharter::hc_add_series(data=ds,type="scatter",#col=cols[1],
                        tooltip = list(headerFormat="",
                          pointFormat=paste("{point.label}","{point.z}")),
                        showInLegend=FALSE)

    hc <- hc %>%
        highcharter::hc_add_series(name = mytypes[1],
		data=line1,type='line',#name=sprintf("regression line %s",type1),
		color = cols[1],enableMouseTracking=FALSE,marker=FALSE) %>%
        highcharter::hc_add_series(name = mytypes[2],
		data=line2,type='line',#name=sprintf("regression line %s",type2),
		color = cols[2],enableMouseTracking=FALSE,marker=FALSE)
    
    hc
}
