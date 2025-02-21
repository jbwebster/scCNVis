#' Create a custom heatmap of copy number data
#'
#' This function takes a SCCNVisObject and outputs a ComplexHeatmap::Heatmap object
#'
#' @param obj A SCCNVisObject made using createPlotObject()
#' @param gr A GenomicRanges:GRanges object of the regions to graph. If NULL, all regions will be plotted. Default = NULL
#' @param add.noise Add minor noise to the data for clustering purposes. Default = TRUE
#' @param color.ramp A color ramp generated using circlize::colorRamp2. If NULL, a red/white/blue color ramp will be used, centered on the median. Default = NULL
#' @param verbose Verbose logging. Default = FALSE
#' @return Returns a ComplexHeatmap object
#' @export
makeSCHeatmap <- function(obj,
                          gr = NULL,
                          add.noise = T,
                          color.map = NULL,
                          verbose = F) {
  .validateObject(obj)

  .verboseLog(verbose, "Running makeSCHeatmap()")

  plot.data <- obj@Matrix
  plot.granges <- obj@GRanges
  if(!is.null(gr)) {
    .verboseLog(verbose, "Subsetting to desired regions")

    overlap <- GenomicRanges::findOverlaps(plot.granges, gr)
    plot.granges <- plot.granges[overlap@from]
    if(length(plot.granges) < 1) {
      stop("Provided GRanges object does not overlap with any regions in the SCCNVisObject")
    }
    plot.granges <- .standardizeGRanges(plot.granges)
    plot.data <- plot.data[,plot.granges$index]
    .verboseLog(verbose, paste0("Genomic bins remaining: ", ncol(plot.data)))
  }

  ###
  # Add noise
  cluster.data <- plot.data
  if (add.noise) {
    .verboseLog(verbose,"Adding noise for clustering")

    noise <- matrix(rnorm(n = length(plot.data), mean = 0, sd = 0.005),
                    nrow = nrow(plot.data),
                    ncol = ncol(plot.data))
    cluster.data <- plot.data + noise
  }

  ###
  # Cluster
  .verboseLog(verbose,"Clustering")
  d <- parallelDist::parallelDist(cluster.data)
  hc_re <- hclust(d, method="ward.D2")

  ###
  # Color scheme
  if(is.null(color.map)) {
    .verboseLog(verbose,"Using default color scheme")

    quantiles <- quantile(plot.data, probs = c(0.01, 0.99))

    cl <- c("Blue","White","Red")
    color.map <- circlize::colorRamp2(
      c(quantiles[1], median(plot.data), quantiles[2]),
      cl
    )
  }

  ####
  #Column annotation and splitting
  .verboseLog(verbose,"Creating column annotation")

  #TODO: May cause odd behavior if granges are unsorted
  n.chr <- length(unique(GenomeInfoDb::seqnames(plot.granges)))
  col.anno.colors <- rep("white",n.chr)
  col.anno.colors[2*c(1:floor(n.chr/2))] <- "lightgray"
  col.annot <- ComplexHeatmap::HeatmapAnnotation(
    foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = col.anno.colors),
                     labels = unique(gsub("chr","",unique(GenomeInfoDb::seqnames(plot.granges)))),
                     labels_gp = grid::gpar(col="black", fontsize=12))
  )
  col.split <- factor(gsub("chr","",GenomeInfoDb::seqnames(plot.granges)),
                      levels=gsub("chr","",unique(GenomeInfoDb::seqnames(plot.granges))))

  ###
  # Plotting
  .verboseLog(verbose, "Plotting")

  ht <- ComplexHeatmap::Heatmap(
    plot.data,
    name = "Heatmap",
    col = color.map,
    cluster_rows = hc_re,
    cluster_columns = FALSE,
    #left_annotation = row.annot,
    top_annotation = col.annot,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    column_split = col.split,
    column_title = NULL,
    cluster_column_slices = FALSE,
    show_column_dend = FALSE,
    use_raster = FALSE
  )

  ht

}


#' Make a line plot showing the average score across different regions
#'
#' @param obj SCCNVis object
#' @param gr A GenomicRanges::GRanges object specifying which reasons to plot. If NULL, all regions will be plotted. Default = NULL
#' @param window.size Window size for moving average calculation
#' @param clusters Character vector specifying clusters for cells. If provided, should be same as the number of cells plotted. Default = NULL
#' @param verbose Verbose logging. Default = FALSE
#' @return Returns a ggplot2::gg object
#' @export
makeLinePlot <- function(obj,
                         gr = NULL,
                         window.size = 5,
                         clusters = NULL,
                         verbose = F) {

  .validateObject(obj)

  .verboseLog(verbose, "Running makeLinePlot()")

  plot.data <- obj@Matrix
  plot.granges <- obj@GRanges
  if(!is.null(gr)) {
    .verboseLog(verbose, "Subsetting to desired regions")

    overlap <- GenomicRanges::findOverlaps(plot.granges, gr)
    plot.granges <- plot.granges[overlap@from]
    if(length(plot.granges) < 1) {
      stop("Provided GRanges object does not overlap with any regions in the SCCNVisObject")
    }
    plot.granges <- .standardizeGRanges(plot.granges)
    plot.data <- plot.data[,plot.granges$index]
    .verboseLog(verbose, paste0("Genomic bins remaining: ", ncol(plot.data)))
  }

  #Calculating column means
  col.means <- data.frame(rn = c(1:ncol(plot.data)), Mean = colMeans(plot.data))

  #Rolling means
  .verboseLog(verbose, paste0("Calculating rolling means, with window size:", window.size))
  moving.average <- .movingAverage(col.means, data.frame(plot.granges), window.size = window.size)

  #Combine rolling means and granges
  mgr <- merge(moving.average,data.frame(plot.granges), by="index", all.y = T)
  #Something this merge causes some NA values for some reason
  mgr[mgr$index %in% moving.average$index,"moving.average"] <- mgr$moving.average

  p <- ggplot2::ggplot(mgr) +
        ggplot2::geom_line(ggplot2::aes(x=index, y=moving.average)) +
        ggplot2::facet_grid(~seqnames, scales = 'free_x', space = 'free_x', switch = 'x') +
        ggplot2::theme_classic() +
        ggplot2::labs(x="",y="") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.margin = ggplot2::unit(0, "lines")) +
        ggplot2::scale_x_continuous(expand = c(0.1, 0.1))

  p

}
