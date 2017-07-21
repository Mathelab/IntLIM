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
	palette = c("#C71585", "#00E5EE")) {

    if ( viewer == TRUE ){
      if (length(palette) == 2) {
        cols <- c(palette)
      }
      else if (length(palette) == 1) {
        cols <- RColorBrewer::brewer.pal(2, palette)
      }
      else {
        stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
      }
    }
    else{
      if(!is.null(palette)){
        cols <- RColorBrewer::brewer.pal(2, palette)
      }
    }
    categ <- c("Genes","Metabolites")

	mygene <- as.data.frame(Biobase::assayDataElement(inputData[["expression"]],'exprs'))
	toplot <- reshape2::melt(mygene)

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
       highcharter::hc_add_series_boxplot(toplot$value,by=toplot$variable,col=cols[1]) %>%
      highcharter::hc_yAxis(title = list(text = "log(expression)",
                            style = list(fontSize = "13px")),
               labels = list(format = "{value}")) %>%
      highcharter::hc_xAxis(labels="") %>%
      highcharter::hc_colors(cols) %>%
      highcharter::hc_tooltip(valueDecimals = 2) %>%
      highcharter::hc_exporting(enabled = TRUE)

	mymetab <- Biobase::assayDataElement(inputData[["metabolite"]],'metabData')
	toplot <- reshape2::melt(mymetab)

        m <- highcharter::highchart(width = 750, height = 750 ) %>%
      highcharter::hc_title(text = "Metabolite Levels",
               style = list(color = '#2E1717',
                            fontWeight = 'bold', fontSize = "20px")) %>%
      highcharter::hc_plotOptions(
        boxplot = boxplotOptions
        ) %>%
       highcharter::hc_add_series_boxplot(toplot$value,by=toplot$variable,col=cols[2]) %>%
      highcharter::hc_yAxis(title = list(text = "log(levels)",
                            style = list(fontSize = "13px")),
               labels = list(format = "{value}")) %>%
      highcharter::hc_xAxis(labels="") %>%
      highcharter::hc_colors(cols[2]) %>%
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
#' @param stype category to color-code by
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
#' PlotPCA(mydata,stype = "PBO_vs_Leukemia")
#' @export
PlotPCA <- function(inputData,viewer=T,stype=NULL,
        palette = c("#C71585", "#00E5EE")) {

    if ( viewer == TRUE ){
      if (length(palette) == 2) {
        cols <- c(palette)
      }
      else if (length(palette) == 1) {
        cols <- RColorBrewer::brewer.pal(2, palette)
      }
      else {
        stop("palette must either be an RColorBrewer palette or a vector of hex colors of size 2")
      }
    }
    else{
      if(!is.null(palette)){
        cols <- RColorBrewer::brewer.pal(2, palette)
      }
    }
    categ <- c("Genes","Metabolites")
    
        if(is.null(stype)) {
		warning("The resulting PCA plot is not color-coded because you did not provide a category in 'stype'")
		mytype <- NULL
        } else if (length(intersect(colnames(Biobase::pData(inputData[["metabolite"]])),stype))!=1) {
		stop(paste0("You provided ",stype, "as your stype variable but id does not exist in your data"))}
        else {
        	mytype <- Biobase::pData(inputData[["metabolite"]])[,stype]
	}
        mygene <- as.data.frame(Biobase::assayDataElement(inputData[["expression"]],'exprs'))
        gpca <- stats::prcomp(t(mygene),center=F,scale=T)
	percvar=round((gpca$sdev)^2 / sum(gpca$sdev^2)*100,2)
        mycols <- as.character(mytype)
	mycols[which(mytype==unique(mytype)[1])] <- cols[1]
	mycols[which(mytype==unique(mytype)[2])] <- cols[2]
 
	if(is.null(mytype)) {
		toplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),color=rep("blue",nrow(gpca$x)))
	} else {
		toplot=data.frame(x=gpca$x[,1],y=gpca$x[,2],z=rownames(gpca$x),label=mytype,color=mycols)
	}
	ds <- list_parse(toplot)

	pg <- highcharter::highchart(width = 350, height = 350 ) %>%
		highcharter::hc_title(text="PCA of genes") %>%
		highcharter::hc_xAxis(title=list(text=paste0("PC1:",round(percvar[1],1),"%"))) %>%
		highcharter::hc_yAxis(title=list(text=paste0("PC2:",round(percvar[2],2),"%"))) %>%
		hc_chart(zoomType = "xy") %>% 
		highcharter::hc_add_series(data=ds,type="scatter",col=cols[1],
			tooltip = list(headerFormat="", 
			  pointFormat=paste("{point.label}","{point.z}")),
			showInLegend=FALSE)
#		dataLabels= list(enabled = TRUE, format = "{point.label}"),

        mymetab <- Biobase::assayDataElement(inputData[["metabolite"]],'metabData')
        mpca <- stats::prcomp(t(mymetab),center=F,scale=T)
	percvar=round((mpca$sdev)^2 / sum(mpca$sdev^2)*100,2)

	if(is.null(mytype)) {
                toplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),color=rep("blue",nrow(mpca$x)))
        } else {
                toplot=data.frame(x=mpca$x[,1],y=mpca$x[,2],z=rownames(mpca$x),label=mytype,color=mycols)
        }
        ds <- list_parse(toplot)

        pm <- highcharter::highchart(width = 350, height = 350 ) %>%
                highcharter::hc_title(text="PCA of metabolites") %>%
                highcharter::hc_xAxis(title=list(text=paste0("PC1:",round(percvar[1],1),"%"))) %>%
                highcharter::hc_yAxis(title=list(text=paste0("PC2:",round(percvar[2],2),"%"))) %>%
                hc_chart(zoomType = "xy") %>%
                highcharter::hc_add_series(data=ds,type="scatter",col=cols[1],
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
    y<-as.numeric(IntLimResults)
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
#' @param metabName string of select metabName
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="IntLim", mustWork=TRUE)
#' csvfile <- file.path(dir, "NCIinput.csv")
#' mydata <- ReadData(csvfile,metabid='id',geneid='id')
#' stype = inputData[["expression"]]$PBO_vs_Leukemia
#' PlotGMPair(inputData,stype,"PRPF8","(p-Hydroxyphenyl)lactic acid")
#' 
#' }
#' @export
PlotGMPair<- function(inputData,stype,geneName,metabName) {
    gene<-Biobase::exprs(inputData[["expression"]])
    sGene<-gene[geneName,]
    
    metab<-Biobase::assayDataElement(inputData[["metabolite"]], 'metabData')
    sMetab<-as.numeric(metab[metabName,])
    
    data<-data.frame('gene'=sGene,'metab'=sMetab,'type'=stype)
    
    
    data<-data[data$type!="",]
    data$type <- factor(data$type)
    b<-stats::glm(data$metab~data$gene+data$type+data$gene:data$type)
    
    
    coefficients<-t(b$coefficients)
    i<-coefficients[,1]
    g<-coefficients[,2]
    c<-coefficients[,3]
    g.c<-coefficients[,4]
    
    u<-as.matrix(levels(data$type))
    type1<-u[1,1]
    type2<-u[2,1]
    
    max<- max(data$gene)
    min<-min(data$gene)
    
    line1<-data.frame(x=c(max,min),y=c(g*max+i,min*g+i))
    line2<-data.frame(x=c(max,min),y=c(g*max+g.c*max+i+c,min*g+i+c))
    
    
    
    hc<-highcharter::highchart(width = 800, height = 700) 
    for(type in u){
        hc<-hc%>% 
            highcharter::hc_add_series_scatter(data$gene[data$type==type],data$metab[data$type==type],name=sprintf("type %s", type),
                                               showInLegend = TRUE) 
        
    } 
    
    hc <- hc %>%
        highcharter::hc_add_series(data=line1,type='line',name=sprintf("regression line %s",type1),color = "#6AB9FF",enableMouseTracking=FALSE,marker=FALSE) %>%
        highcharter::hc_add_series(data=line2,type='line',name=sprintf("regression line %s",type2),color = "#474544",enableMouseTracking=FALSE,marker=FALSE)
    
    hc
}
