
# SCCNVis: Visualizing single-cell copy number data

SCCNVis is a small project for visualizing copy number data from single-cell experiments. The project mostly serves as a wrapper for ComplexHeatmap and ggplot2.

Multiple tools exist for calling copy number variants using single cell data. Most of those tools include default plotting methods which work well in most cases, especially since many users want to simply use the copy number calls as a way to confirm the cancer status of individual cells. However, in situations where deeper exploration of single-cell copy number results are desired, the basic plotting methods of those tools are insufficient. This package contains helper functions as a starting point to create additional visualizations to explore copy number data.

## Installation
`devtools::install_github("jbwebster/SCCNVis")`

Some dependencies may fail to install because they must be installed through BiocManager rather than standard CRAN. If that is the case, you can try to install those dependencies by doing the following:

`
if (!require("BiocManager", quietly = TRUE) )
	install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
BiocManager:;install("GenomicRanges")
`

## Tutorial

A high-level overview of how the package works is outlined below. More in-depth examples can be found HERE.

SCCNVis works with SCCNVis objects. This provides a common format that can be created using the outputs from any of the common CNV callers. To create a SCCNVis object, you need 1) a matrix (rows = cells, columns = regions) with copy number values, 2) a list of cell names (should be the same length as nrows(matrix) and 3) a GenomicRanges::GRanges object describing the columns of the matrix. For example, if the input matrix describes 1000 cells and a binned genome resulting in 3000 bins, then the cells list should have 1000 cell names and the GRanges object should describe 3000 genomic regions. The package includes some example data, which I use below

`
library(SCCNVis)
input.matrix <- example.matrix
cell.names <- example.cells
granges <- example.granges
meta <- example.meta
obj <- createPlotObject(input.matrix, cell.names, granges, example.meta)
`

The resulting `obj` output can then be passed into any of the available plotting functions. For example:

`
heatmap.result <- makeSCHeatmap(obj)
`

The output of the plotting functions is a list with two items. The first is the Plot (which can be accessed as heatmap.result$Plot) and the second is the data that was used to make the plot (accessed as heatmap.result$PlotData) in case you wish to create your own version of the plot or inspect the data more thoroughly. If you wish to simply export the plot, it can be done using the wrapper function:

`
saveCustomPlot(heatmap.result, "~/path/to/heatmap.png")
`


