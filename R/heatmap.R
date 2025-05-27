#' Create a custom heatmap of copy number data
#'
#' This function takes a SCCNVisObject and outputs a ComplexHeatmap::Heatmap object
#'
#' @param obj A SCCNVisObject made using createPlotObject()
#' @param gr A GenomicRanges:GRanges object of the regions to graph. If NULL, all regions will be plotted. Default = NULL
#' @param annotations Column names from obj metadata to use for annotating cells on the left of the plot. Default = NULL
#' @param show.annotation.name Show annotation names. Does nothing if annotations is NULL. Default == NULL.
#' @param add.noise Add minor noise to the data for clustering purposes. The visualization will not include the noise. Default = TRUE
#' @param color.ramp A color ramp generated using circlize::colorRamp2. If NULL, a red/white/blue color ramp will be used, centered on the median. Default = NULL
#' @param legend.title Title for the legend describing the heatmap. Often something like "log2(copy ratio)" or some other metric. Default = "Value"
#' @param verbose Verbose logging. Default = FALSE
#' @return Returns a list containing a "Plot" (ComplexHeatmap) and "PlotData" (matrix)
#' @export
makeSCHeatmap <- function(obj,
                          gr = NULL,
                          annotations = NULL,
                          annotation.colors = NULL,
                          show.annotation.name = FALSE,
                          add.noise = T,
                          color.map = NULL,
                          legend.title = "Value",
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


  #Annotation
  if(!is.null(annotations)) {
    if (sum(annotations %in% colnames(obj@Meta)) == length(annotations)) {
      annotation_list <- list()
      for (annot in annotations) {
        annotation_list[[annot]] <- obj@Meta[,annot]
      }
      if(is.null(annotation.colors)) {
        args <- c(annotation_list, list(
          which = "row",
          show_annotation_name = show.annotation.name,
        ))
      } else {
        args <- c(annotation_list, list(
          which = "row",
          show_annotation_name = show.annotation.name,
          col = annotation.colors
        ))
      }

      row.annot <- do.call(ComplexHeatmap::HeatmapAnnotation, args)

    }
  } else {
    row.annot <- NULL
  }


  ###
  # Plotting
  .verboseLog(verbose, "Plotting")

  ht <- ComplexHeatmap::Heatmap(
    plot.data,
    name = legend.title,
    col = color.map,
    cluster_rows = hc_re,
    cluster_columns = FALSE,
    left_annotation = row.annot,
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

  return(list("Plot" = ht, "PlotData" = list("PlotData" = plot.data, "Clustering" = hc_re)))

}

