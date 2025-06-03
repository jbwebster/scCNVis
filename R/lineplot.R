

#' Make a line plot showing values across different regions
#'
#' @param obj scCNVis object
#' @param gr A GenomicRanges::GRanges object specifying which regions to plot. If NULL, all regions will be plotted. Default = NULL
#' @param window.size Window size for moving average calculation
#' @param group Name of metadata column to group lines by. Default = NULL
#' @param group.colors Named character vector. Names should be the values in the column specified by 'group', values are the colors. Default = NULL
#' @param group.label The label to be used for the groups. Default = ""
#' @param method Method for aggregating values across a group. Should be "mean", "median", "sd" or "range". Default = "mean"
#' @param plot.average.of.groups Logical. If TRUE, plot each cluster individually and an additional line representing the average of each cluster. Default = FALSE
#' @param avg.color Color. Used for the color of the cluster average, if plot.average.of.clusters == T. Default = 'red'
#' @param show.x.ticks Logical. If TRUE, show tick marks on X axis. Only recommended when showing smaller regions. Default = FALSE
#' @param ylab Character. Label of y axis. Default = ""
#' @param remove.chr.prefix Remove 'chr' prefix of chromosome names. Default = TRUE
#' @param verbose Verbose logging. Default = FALSE
#' @return Returns a list containing "Plot" (ggplot) and "PlotData" (data.frame)
#' @export
makeLinePlot <- function(obj,
                         gr = NULL,
                         window.size = 5,
                         group = NULL,
                         group.colors = NULL,
                         group.label = "",
                         method = "mean",
                         plot.average.of.groups = FALSE,
                         avg.color = 'red',
                         show.x.ticks = FALSE,
                         ylab = "",
                         remove.chr.prefix = TRUE,
                         verbose = F) {

  .validateObject(obj)

  .verboseLog(verbose, "Running makeLinePlot()")

  plot.data <- obj@Matrix
  plot.granges <- obj@GRanges
  if(!is.null(gr)) {
    .verboseLog(verbose, "Subsetting to desired regions")

    overlap <- GenomicRanges::findOverlaps(plot.granges, gr, minoverlap = 10)
    plot.granges <- plot.granges[overlap@from]
    if(length(plot.granges) < 1) {
      stop("Provided GRanges object does not overlap with any regions in the scCNVisObject")
    }
    plot.granges <- .standardizeGRanges(plot.granges)
    plot.data <- plot.data[,plot.granges$index]
    .verboseLog(verbose, paste0("Genomic bins remaining: ", ncol(plot.data)))
  }

  plot.groups <- NULL
  if (!is.null(group)) {
    if (!(group %in% colnames(obj@Meta))) {
      stop("The value in 'group' is not a column in obj@Meta")
    } else {
      plot.groups <- obj@Meta[,group]
      names(plot.groups) <- rownames(obj@Meta)
    }
  } else {
    plot.groups <- rep("All cells", nrow(plot.data))
    names(plot.groups) <- rownames(plot.data)
  }

  res.mgr <- NULL
  for (curr.group in unique(plot.groups)) {
    cells <- plot.groups[plot.groups==curr.group]
    .verboseLog(verbose, paste0("Calculating moving average for group: ", curr.group, " using window size: ", window.size))
    curr.data <- plot.data[names(cells),]
    if (method == "mean") {
      col.means <- data.frame(rn = plot.granges$index, Mean = colMeans(curr.data))
    } else if (method == "median") {
      col.means <- data.frame(rn = plot.granges$index, Mean = apply(curr.data,2,median))
    } else if (method == "range") {
      col.means <- data.frame(rn = plot.granges$index, Mean = apply(curr.data,2,max) - apply(curr.data,2,min))
    } else if (method == "sd") {
      col.means <- data.frame(rn = plot.granges$index, Mean = apply(curr.data,2,sd))
    } else {
      stop("method must be 'mean', 'median', 'sd', or 'range'")
    }

    moving.average <- .movingAverage(col.means, data.frame(plot.granges), window.size = window.size)
    mgr <- merge(moving.average,data.frame(plot.granges), by="index", all.y = TRUE)
    mgr[mgr$index %in% moving.average$index,"moving.average"] <- mgr$moving.average
    mgr$group <- rep(curr.group,nrow(mgr))
    if (is.null(res.mgr)) {
      res.mgr <- mgr
    } else {
      res.mgr <- rbind(res.mgr, mgr)
    }
  }
  res.mgr$group <- factor(res.mgr$group)
  res.mgr$middle <- res.mgr$start + ((res.mgr$end - res.mgr$start) / 2)
  if(remove.chr.prefix) {
    res.mgr$seqnames <- gsub("chr","",res.mgr$seqnames)
    res.mgr$seqnames <- factor(res.mgr$seqnames,
                               levels = unique(res.mgr$seqnames),
                               ordered = TRUE)

  }


  if(plot.average.of.groups) {
    grouped <- res.mgr %>%
      dplyr::group_by(index, seqnames, start, end, middle, width, strand) %>%
      dplyr::summarise(moving.average = mean(moving.average), group = "Grouped.Mean")
    res.mgr <- rbind(res.mgr,grouped)
    res.mgr$group <- factor(res.mgr$group)

    p <- ggplot2::ggplot(data=res.mgr) +
      ggplot2::geom_line(data=res.mgr[res.mgr$group != "Grouped.Mean",],
                         ggplot2::aes(x=middle,y=moving.average,group=group),
                         alpha=0.25,size=0.5,color="black") +
      ggplot2::geom_line(data=res.mgr[res.mgr$group == "Grouped.Mean",],
                         ggplot2::aes(x=middle,y=moving.average,group=group),
                         alpha=1,size=2,color=avg.color) +
      ggplot2::facet_grid(~seqnames, scales = "free_x", space = "free_x") +
      ggplot2::theme_classic() +
      ggplot2::labs(x="",y=ylab) +
      ggplot2::theme(panel.spacing = ggplot2::unit(0,"lines"))
  } else {

    p <- ggplot2::ggplot(res.mgr) +
      ggplot2::geom_line(ggplot2::aes(x=middle, y=moving.average, color=group)) +
      ggplot2::facet_grid(~seqnames, scales = 'free_x', space = 'free_x') +
      ggplot2::theme_classic() +
      ggplot2::labs(x="",y=ylab,color=group.label) +
      ggplot2::theme(panel.spacing = ggplot2::unit(0, "lines"))

    if (!is.null(group.colors)) {
      if (length(group.colors) == length(unique(res.mgr$group))) {
        p + ggplot2::scale_color_manual(values = group.colors)
      } else {
        message(paste0("group.colors length is different than the number of groups present (", length(unique(res.mgr$group)), ")"))
        stop()
      }
    }
  }

  if (show.x.ticks) {
    p <- p + ggplot2::scale_x_continuous(expand = c(0.1,0.1),
                                    labels = scales::unit_format(unit = "Mb",  scale = 1e-6))
  } else {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::scale_x_continuous(expand = c(0.1, 0.1))

  }




  return(list("Plot" = p, "PlotData" = res.mgr))

}
