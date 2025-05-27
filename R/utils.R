

#' Convert AtaCNV column names into a GRanges object
#'
#' This functions converts the column names provided by AtaCNV output into
#' a GenomicRanges::GRanges object for use in creating a SCCNVis object.
#' Adapted from Aelita-Stone/AtaCNV
#'
#' @param bin.names Character vector of column names from AtaCNV output.
#' @return GenomicRanges::GRanges object
#' @export
ataCNVToGRanges <- function(bin.names) {
  temp <- strsplit(bin.names, split = "_")
  temp2 <- unlist(temp)
  dim(temp2) <- c(3,length(temp2)/3)
  bins <- data.frame(chr=temp2[1,],
                     start=as.numeric(temp2[2,]),
                     end=as.numeric(temp2[3,]),
                     strand=rep("+",length(temp2[1,])))
  g <- GenomicRanges::makeGRangesFromDataFrame(bins,
                                               seqnames.field = "chr",
                                               start.field = "start",
                                               end.field = "end",
                                               strand.field = "strand")
  g <- .standardizeGRanges(g)
  g
}


#' Add a column to the metadata held in the SCCNVis object
#'
#' @param obj SCCNVis object to add metadata to
#' @param meta Vector of values to be added to the metadata
#' @param colname Name of new column to add to metadata
#' @param cells Cells to add metadata to. If NULL, it is assumed that meta contains one value for each cell and is in the same order as the object's current metadata. If provided, it should be in the same order as meta. Default = NULL
#' @return SCCNVis object with updated metadata
#' @export
addMetaData <- function(obj,meta,colname,cells=NULL) {
  if (class(obj) != "SCCNVisObject") {
    stop("obj is not a SCCNVisObject")
  }

  if (is.null(meta)) {
    stop("meta value(s) cannot be missing")
  }

  if (is.null(colname)) {
    stop("colname can not be NULL")
  }

  if(sum(cells %in% obj@Meta$CellNames) < length(cells)) {
    stop("Not all provided cell names are in obj@Meta$CellNames")
  }

  if(sum(duplicated(cells)) > 0) {
    stop("Duplicate cell names provided")
  }

  if (is.null(cells)) {
    if (length(meta) != nrow(obj@Meta)) {
      stop("Length of meta is expected to equal nrow(obj@Meta) when cells == NULL")
    }
    df <- data.frame(CellNames = obj@Meta$CellNames,
                     X = meta)
    colnames(df) <- c("CellNames", colname)
    merged <- merge(obj@Meta, df, by="CellNames")
    rownames(merged) <- merged$CellNames
    merged <- merged[rownames(obj@Meta),]
    obj@Meta <- merged
    return(obj)
  } else {
    names(meta) <- cells

    df <- data.frame(CellNames = obj@Meta$CellNames,
                     X = rep(NA, nrow(obj@Meta)))
    keys <- df[df$CellNames %in% cells,"CellNames"]
    df[df$CellNames %in% cells,"X"] <- meta[keys]
    colnames(df) <- c("CellNames",colname)
    merged <- merge(obj@Meta, df, by="CellNames")
    rownames(merged) <- merged$CellNames
    merged <- merged[rownames(obj@Meta),]
    obj@Meta <- merged
    return(obj)

  }
}


#' Create a SCCNVis object for plotting
#'
#' This function creates a basic object containing the information needed
#' to plot copy number data from a single-cell experiment. It expects
#' an input matrix where each row represents a cell and each column represents
#' a region of the genome.
#'
#' @param input.matrix An input matrix of class 'matrix'
#' @param cell.names A character vector containing unique cell names. Should be equal to the number of rows in input.matrix
#' @param granges A GenomicRanges GRanges object, where each range in the object corresponds to a column in the input.matrix
#' @param meta Data.frame of meta data. Rownames should be cell.names. Default = NULL
#' @return A list object containing 2 elements ("Matrix" and "GRanges")
#' @export
createPlotObject <- function(input.matrix, cell.names, granges, meta = NULL) {
  .validObjectInputs(input.matrix, cell.names, granges, meta)

  ###
  #Handle input.matrix
  input.matrix <- as.matrix(input.matrix)

  ###
  #Doesn't modify cell.names unless necessary
  #Warning is given in .validObjectInputs
  cell.names <- make.unique(cell.names)

  rownames(input.matrix) <- cell.names

  ###
  #Handle GRanges
  granges <- .standardizeGRanges(granges)
  granges$index <- c(1:length(granges))

  ###
  #Handle metadata
  meta <- .standardizeMeta(meta, cell.names)


  out <- new("SCCNVisObject",
             Matrix = input.matrix,
             GRanges = granges,
             Meta = meta)
}



#' Helper function for saving a generated plot
#'
#' A very basic helper function for saving plots with default settings.
#' If more customized settings are necessary, it is recommended to see methods
#' associated with the used plotting packages, as this is a simplified wrapper
#' for those methods. In the case of Heatmaps, see
#' ComplexHeatmap::draw(). For all other plots, see ggplot2::ggsave().
#'
#' @param p Plot object
#' @param filename File path and filename for output
#' @param width Width of plot, in inches
#' @param height Height of plot, in inches
#' @param format Format of output file. Should be "png" or "pdf"
#' @return None
#' @export
saveCustomPlot <- function(p, filename, width = 14, height = 7, format = "png") {
  if (format != "png" & format != "pdf") {
    stop("format must == 'png' or 'pdf'")
  }

  if (is.null(filename)) {
    stop("filename must not be NULL")
  }

  #Non-heatmaps
  if(length(p)>1) {
    if (format == "png") {
      ggplot2::ggsave(filename, width = width, height = height, units = "in", dpi = 700)
    }
    if (format == "pdf") {
      ggplot2::ggsave(filename, width = width, height = height, units = "in", dpi = 700)
    }
  }
  else if(class(p) == "Heatmap") {
    if (format == "png") {
      png(filename, width = width, height = height, units = "in", res = 1000)
      ComplexHeatmap::draw(p)
      dev.off()
    }
    if (format == "pdf") {
      pdf(filename, width = width, height = height)
      ComplexHeatmap::draw(p)
      dev.off()
    }
  }
}


.movingAverage <- function(input, df.gr, window.size = 5) {
  out <- data.frame("index" = 0, "moving.average" = 0)
  for(seqname in unique(df.gr$seqnames)) {
    curr.indices <- df.gr[df.gr$seqnames==seqname,"index"]
    curr.mat <- input[input$rn %in% curr.indices,]
    ma <- stats::filter(curr.mat, rep(1 / window.size, window.size), sides = 2)
    out <- rbind(out, data.frame("index" = curr.indices, "moving.average" = ma[,2]))
  }
  out[2:nrow(out),]
}


.verboseLog <- function(v, msg) {
  if(v) {
    message(msg)
  }
}

