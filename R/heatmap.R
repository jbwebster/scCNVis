#' Create a custom heatmap of copy number data
#'
#' This function takes a SCCNVisObject and outputs a ComplexHeatmap::Heatmap object
#'
#' @param obj A SCCNVisObject made using createPlotObject()
#' @param gr A GenomicRanges:GRanges object of the regions to graph. If NULL, all regions will be plotted. Default = NULL
#' @param annotations Column names from obj metadata to use for annotating cells on the left of the plot. Default = NULL
#' @param annotation.colors
#' @param group Name of metadata column. Signal will be averaged across all cells in each group. Default = NULL
#' @param secondary.group Name of metadata column. Can be used to add a secondary annotation describing the composition of 'group'. Default = NULL
#' @param show.annotation.name Show annotation names. Does nothing if annotations is NULL. Default == NULL
#' @param add.noise Add minor noise to the data for clustering purposes. The visualization will not include the noise. Default = TRUE
#' @param heatmap.colors Named vector of numeric values. Names should be colors, values correspond to values in the heatmap. If NULL, a red/white/blue color ramp will be used, centered on the median. Default = NULL
#' @param legend.title Title for the legend describing the heatmap. Often something like "log2(copy ratio)" or some other metric. Default = "Value"
#' @param verbose Verbose logging. Default = FALSE
#' @return Returns a list containing a "Plot" (ComplexHeatmap) and "PlotData" (matrix)
#' @export
makeSCHeatmap <- function(obj,
                          gr = NULL,
                          annotations = NULL,
                          annotation.colors = NULL,
                          group = NULL,
                          secondary.group = NULL,
                          show.annotation.name = FALSE,
                          add.noise = T,
                          heatmap.colors = NULL,
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
  # Calculate group averages, if desired
  if(!is.null(group)) {
    if(group %in% colnames(obj@Meta)) {
      groups <- obj@Meta[,group]
      names(groups) <- rownames(obj@Meta)
      groups <- groups[!is.na(groups)]
      k <- length(unique(groups))
      grouped.mat <- matrix(nrow = k, ncol = ncol(plot.data))
      rownames(grouped.mat) <- unique(groups)[order(unique(groups))]
      for (kc in order(unique(groups))) {
        grouped.mat[kc,] <- colMeans(plot.data[groups==kc,], na.rm = TRUE)
      }
      d <- parallelDist::parallelDist(grouped.mat)
      hc_re <- hclust(d, method="ward.D2")
      new.kclusters <- cutree(hc_re, k)
      plot.data <- grouped.mat
    } else {
      stop("group must be a column name in obj@Meta")
    }
  } else {
    k <- NULL
  }


  ###
  # Color scheme
  if(is.null(heatmap.colors)) {
    .verboseLog(verbose,"Using default color scheme")

    quantiles <- quantile(plot.data, probs = c(0.01, 0.99))

    cl <- c("Blue","White","Red")
    color.map <- circlize::colorRamp2(
      c(round(quantiles[1],2), round(median(plot.data),2), round(quantiles[2],2)),
      cl
    )
    heatmap.colors <- c(round(quantiles[1],2), round(median(plot.data),2), round(quantiles[2],2))
    names(heatmap.colors) <- cl
  } else {
    color.map <- circlize::colorRamp2(
      heatmap.colors,
      names(heatmap.colors)
    )
  }

  #Make default legend using custom color map
  default_legend <- ComplexHeatmap::Legend(
    title = legend.title,
    at = heatmap.colors,
    col_fun = color.map
  )


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
  .verboseLog(verbose,"Creating left-side annotation, if any")
  if(!is.null(annotations) & is.null(k)) {
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
  } else if(!is.null(k)) {
    max.color <- 2 + k
    row.annot <- ComplexHeatmap::HeatmapAnnotation(
      Subclone = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 3:max.color),
                            labels = rownames(plot.data)),
      which = "row"
    )
  } else {
    row.annot <- NULL
  }

  #Right side annotation, only used if k is not NULL
  if(!is.null(k)) {
    ncells <- summary(factor(groups))
    ncells <- ncells[rownames(plot.data)]
    if(is.null(secondary.group)) {
      right.row.annot <- ComplexHeatmap::HeatmapAnnotation(
        nCells = ComplexHeatmap::anno_barplot(ncells,
                                            baseline = 0,
                                            add_numbers = TRUE,
                                            width = grid::unit(2, "cm"),
                                            gp = grid::gpar(fill = 2, fontsize = 12),
                                            axis_param = list(
                                              gp = grid::gpar(fontsize = 12),
                                              side = "top")
                                            ),
        which = "row",
        annotation_name_rot = 45,
        annotation_label = c("# Cells"),
        annotation_name_gp = grid::gpar(fontsize = 12),
        annotation_name_side = "top"
        )
      bar_legend <- NULL
    } else if (secondary.group %in% colnames(obj@Meta)) {
      secondary.groups <- unique(obj@Meta[,secondary.group])
      prop.matrix <- matrix(rep(0, k* length(secondary.groups)),
                            nrow = k,
                            ncol = length(secondary.groups)
      )
      colnames(prop.matrix) <- secondary.groups
      rownames(prop.matrix) <- rownames(plot.data)
      for (kc in rownames(prop.matrix)) {
        curr.meta <- obj@Meta[obj@Meta[,group] == kc,]
        total <- nrow(curr.meta)
        for (curr.secondary in secondary.groups) {
          curr.x <- curr.meta[curr.meta[,secondary.group] == curr.secondary,]
          curr.n <- nrow(curr.x)
          prop.matrix[kc,curr.secondary] <- curr.n / total
        }
      }
      cl <- length(secondary.groups) + 2
      right.row.annot <- ComplexHeatmap::HeatmapAnnotation(
        nCells = ComplexHeatmap::anno_barplot(ncells,
                                              baseline = 0,
                                              add_numbers = TRUE,
                                              width = grid::unit(2, "cm"),
                                              gp = grid::gpar(fill = 2, fontsize = 12),
                                              axis_param = list(
                                                gp = grid::gpar(fontsize = 12),
                                                side = "top")
        ),
        SecondaryGroup = ComplexHeatmap::anno_barplot(prop.matrix,
                                                      width = grid::unit(2,"cm"),
                                                      gp = grid::gpar(fill = 3:cl, fontsize=12),
                                                      axis_param = list(
                                                        grid::gpar(fontsize = 12),
                                                        side = "top"
                                                      )),
        which = "row",
        annotation_name_rot = 45,
        annotation_label = c("# Cells", secondary.group),
        annotation_name_gp = grid::gpar(fontsize = 12),
        annotation_name_side = "top"
      )

      bar_legend <- ComplexHeatmap::Legend(
        labels = colnames(prop.matrix),
        title = secondary.group,
        legend_gp = grid::gpar(fill=3:cl)
      )
    }

  } else {
    right.row.annot <- NULL
    bar_legend <- NULL
  }



  ###
  # Plotting
  .verboseLog(verbose, "Plotting")

  if(is.null(k)) {
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
  } else {
    ht <- ComplexHeatmap::Heatmap(
      plot.data,
      name = legend.title,
      col = color.map,
      cluster_rows = FALSE,
      row_split = new.kclusters,
      border = TRUE,
      row_title = NULL,
      cluster_columns = FALSE,
      left_annotation = row.annot,
      right_annotation = right.row.annot,
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
  }
  if(is.null(bar_legend)) {
    legend.list <- list(default_legend)
  } else {
    legend.list <- list(default_legend, bar_legend)
  }


  return(list("Plot" = ht, "PlotData" = list("PlotData" = plot.data, "Clustering" = hc_re, "Legend" = legend.list)))

}

