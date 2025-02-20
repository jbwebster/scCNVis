
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
#' @return A list object containing 2 elements ("Matrix" and "GRanges")
#' @export
createPlotObject <- function(input.matrix, cell.names, granges) {
  .validObjectInputs(input.matrix, cell.names, granges)

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

  out <- new("SCCNVisObject",
             Matrix = input.matrix,
             GRanges = granges)
}
