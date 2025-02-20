

.validObjectInputs <- function(input.matrix, cell.names, granges) {
  ####
  #Validate input.matrix
  if(is.null(input.matrix)) {
    stop("Input matrix is null")
  }
  if(sum(class(input.matrix) %in% "matrix") == 0) {
    stop("Input matrix must be of class 'matrix'")
  }

  ####
  #Validate cell.names
  if(is.null(cell.names)) {
    stop("cell.names is null")
  }

  if(class(cell.names) != "character") {
    stop("cell.names must be of class 'character'")
  }

  if(length(cell.names) != nrow(input.matrix)) {
    stop("cell.names must be the same length as the number of rows in input.matrix")
  }

  if(length(unique(cell.names)) != length(cell.names) ) {
    message("Warning: cell.names are not unique. They will be modified to make them unique.")
  }

  ####
  #Validate granges
  if(is.null(granges)) {
    stop("granges is null")
  }

  if(class(granges) != "GRanges") {
    stop("granges must be of class 'GRanges' from the GenomicRanges package")
  }

  if(length(granges) != ncol(input.matrix)) {
    stop("The length of granges should be the same as the number of columns in input.matrix")
  }
}


.validateObject <- function(obj) {
  if(is.null(obj)) {
    stop("Input object is null")
  }

  if(class(obj) != "SCCNVisObject") {
    stop("Input object is not of class SCCNVisObject")
  }
}


.standardizeGRanges <- function(granges) {
  clean.chrm.names <- gsub("[Cc]hr[m]?","",GenomeInfoDb::seqnames(granges))
  clean.chrm.names <- paste0("chr",clean.chrm.names)
  out <- GenomicRanges::GRanges(seqnames = clean.chrm.names,
                 ranges = IRanges::ranges(granges, use.mcols=T),
                 strand = BiocGenerics::strand(granges))

  out
}
