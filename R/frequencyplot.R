
#' Create a frequency plot
#'
#' This function takes a scCNVisObject and outputs a ggplot2 object
#'
#' @param obj A scCNVisObject made using createPlotObject()
#' @param gr A GenomicRanges:GRanges object of the regions to graph. If NULL, all regions will be plotted. Default = NULL
#' @param group Metadata column. Average signal will be calculated per group and then called as loss/neutral/gain. Default = "Sample"
#' @param gain.color Color to represent gain frequency. Default = 'cornflowerblue'
#' @param loss.color Color to represent loss frequency. Default = 'deeppink'
#' @param secondary.group Secondary group for plot groupings. Each value in group, should overlap 1 and only 1 secondary grouping. Default = NULL
#' @param data.func A function for defining thresholds of loss(0)/neutral(1)/gain(2) and otherwise manipulating plot data. If not provided, mean +/- sd will be used for thresholding. See FrequencyPlot example on Github for examples. Default = NULL
#' @param remove.chr.prefix Remove 'chr' prefix of chromosome names. Default = TRUE
#' @param verbose Verbose logging. Default = FALSE
#' @return Returns a list containing a "Plot" (ComplexHeatmap) and "PlotData" (matrix)
#' @export
makeFrequencyPlot <- function(obj,
                          gr = NULL,
                          group = "Sample",
                          gain.color = 'cornflowerblue',
                          loss.color = 'deeppink',
                          secondary.group = NULL,
                          data.func = NULL,
                          remove.chr.prefix = TRUE,
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
      stop("Provided GRanges object does not overlap with any regions in the scCNVisObject")
    }
    plot.granges <- .standardizeGRanges(plot.granges)
    plot.data <- plot.data[,plot.granges$index]
    .verboseLog(verbose, paste0("Genomic bins remaining: ", ncol(plot.data)))
  }

  ##
  # Calculate group means
  if(group %in% colnames(obj@Meta)) {
    groups <- obj@Meta[,group]
    names(groups) <- rownames(obj@Meta)
    groups <- groups[!is.na(groups)]
    k <- length(unique(groups))
    grouped.mat <- matrix(nrow = ncol(plot.data), ncol = k)
    colnames(grouped.mat) <- unique(groups)[order(unique(groups))]
    rownames(grouped.mat) <- colnames(plot.data)
    #Use an index instead of the group values, in case the group values are numeric
    i <- 1
    for (kc in unique(groups)[order(unique(groups))]) {
      grouped.mat[,i] <- colMeans(plot.data[groups==kc,], na.rm = TRUE)
      i <- i + 1
    }
  } else {
    stop("group must be a column in obj@Meta")
  }
  grouped.mat <- data.frame(grouped.mat)

  ##
  # Add regional data and convert to long format
  start.cols <- ncol(grouped.mat)
  gr.df <- data.frame(plot.granges)
  grouped.mat$chrm <- gr.df$seqnames
  grouped.mat$middle <- gr.df$start + (0.5 & gr.df$width)
  grouped.mat$pos.index <- gr.df$index
  plot.data.long <- tidyr::pivot_longer(data = grouped.mat,
                                        cols = all_of(c(1:start.cols)))

  ##
  # Perform thresholding on groups
  if(is.null(data.func)) {
    data.func <- function(data) {
      sd.value <- sd(data$value, na.rm = TRUE)
      mean.value <- mean(data$value, na.rm = TRUE)
      data[data$value > mean.value + sd.value,"value"] <- 2
      data[data$value < mean.value - sd.value,"value"] <- 0
      data[data$value != 2 & data$value != 0,"value"] <- 1
      return(data)
    }
  }

  plot.data.long <- data.func(plot.data.long)
  valid.cols <- c("chrm","middle","pos.index","name","value")
  if (!all(valid.cols %in% colnames(plot.data.long))) {
    stop("data.func output should not remove or rename the columns of the dataframe it is provided")
  }
  if (!is.null(secondary.group)) {
    if(!(secondary.group %in% colnames(obj@Meta))) {
      stop("If provided, secondary.group must be a column in obj@Meta")
    }
    tmp.df <- data.frame(name = obj@Meta[,group],
                         secondary = obj@Meta[,secondary.group])
    tmp.df <- tmp.df[!duplicated(tmp.df),]
    nm <- tmp.df$name
    tmp.df$name <- make.names(nm) #Ensure they get validated the same was as plot.data.long was
    merged <- merge(plot.data.long,tmp.df,by="name")
    
    grouped <- merged %>%
      dplyr::group_by(chrm,middle,pos.index,secondary) %>%
      dplyr::summarise(n_groups = length(unique(name)),
                       n_gain = sum(value == 2),
                       n_loss = sum(value == 0))
    
    #grouped <- plot.data.long %>%
    #  dplyr::group_by(chrm,middle,pos.index,.data[[secondary.group]]) %>%
    #  dplyr::summarise(n_groups = length(unique(name)),
    #                   n_gain = sum(value == 2),
    #                   n_loss = sum(value == 0))
  } else {
    grouped <- plot.data.long %>%
      dplyr::group_by(chrm,middle,pos.index) %>%
      dplyr::summarise(n_groups = length(unique(name)),
                       n_gain = sum(value == 2),
                       n_loss = sum(value == 0))
  }

  grouped$FreqGain <- grouped$n_gain / grouped$n_groups
  grouped$FreqLoss <- -1 * (grouped$n_loss / grouped$n_groups)

  if(remove.chr.prefix) {
    grouped$chrm <- gsub("chr","",grouped$chrm)
    grouped$chrm <- factor(grouped$chrm,
                           levels = unique(grouped$chrm),
                           ordered = TRUE)
  }

  ##
  # Plot
  if(is.null(secondary.group)) {
    p <- ggplot2::ggplot(grouped) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=0,ymax=FreqGain,x=middle),
                           fill=gain.color) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=FreqLoss,ymax=0,x=middle),
                           fill=loss.color) +
      ggplot2::facet_grid( ~ chrm, scales = 'free_x', space = 'free_x') +
      ggplot2::labs(x='Chromosome',y='Frequency of Gain/Loss') +
      ggplot2::theme_classic() +
      ggplot2::theme(panel.spacing = ggplot2::unit(0, "lines"),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
  } else {
    p <- ggplot2::ggplot(grouped) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=0,ymax=FreqGain,x=middle),
                           fill=gain.color) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=FreqLoss,ymax=0,x=middle),
                           fill=loss.color) +
      ggplot2::geom_hline(yintercept=c(-1,1),color='black') +
      ggplot2::scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
      ggplot2::facet_grid(secondary ~ chrm, scales = 'free_x', space = 'free_x') +
      ggplot2::labs(x='Chromosome',y='Frequency of Gain/Loss') +
      ggplot2::theme_classic() +
      ggplot2::theme(panel.spacing.x = ggplot2::unit(0, "lines"),
                     panel.spacing.y = ggplot2::unit(0.5,"lines"),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
  }

  return(list("Plot" = p,
              "PlotData" = grouped))


}
