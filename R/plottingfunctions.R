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
#' csvfile <- file.path(dir, "test.csv")
#' mydata <- ReadData(csvfile,metabid='BIOCHEMICAL',geneid='X')
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
      highcharter::hc_add_series_boxplot(toplot$value,by=toplot$variable,
		col=cols[1], showInLegend = FALSE,outliers=FALSE) %>%
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
      highcharter::hc_add_series_boxplot(toplot$value,by=toplot$variable,
		showInLegend = FALSE,outliers=FALSE) %>%
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


return(p)
}
