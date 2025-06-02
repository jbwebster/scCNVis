
# scCNVis: Single-Cell Copy Number Visualization

scCNVis is a visualization tool for single-cell copy number data.

Multiple tools exist for calling copy number variants using single cell data.
Most of those tools include default plotting methods which work well in most
cases, especially since many users want to simply use the copy number calls as
a way to confirm the cancer status of individual cells. However, in situations
where a deeper exploration of single-cell copy number results is desired,
the basic plotting methods of those tools are insufficient. This package
contains helper functions to support exploration of copy number data and
gain novel biological insights.

Developed by Jace Webster while working in the lab of David Quigley at UCSF.

View the [full documenation here](https://jbwebster.github.io/scCNVis/)

## Installation
```
devtools::install_github("jbwebster/scCNVis")
```

Some dependencies may fail to install because they must be installed through
BiocManager rather than standard CRAN. If that is the case, you can try to
install those dependencies by doing the following:

``` 
if (!require("BiocManager", quietly = TRUE) )
	install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
BiocManager::install("GenomicRanges")
```

## Tutorial

A high-level overview of how the package works is outlined below.

scCNVis works with scCNVis objects. This provides a common format that can be
created using the outputs from any of the common CNV callers. To create a 
scCNVis object, you need:

1. A matrix (rows = cells, columns = regions) with copy number values
2. A list of cell names (should be the same length as nrows(matrix))
3. A GenomicRanges::GRanges object describing the columns of the matrix
4. A metadata dataframe that minimally has rownames as cell names and a column named "Sample"

You can start with the provided example data to get a sense of what these
inputs look like:

```
library(scCNVis)
input.matrix <- example.matrix
cell.names <- example.cells
granges <- example.granges
meta <- example.meta
obj <- scCNVis::createPlotObject(input.matrix, cell.names, granges, meta)
```

The resulting `obj` output can then be passed into any of the available plotting
functions. For example:

```
heatmap.result <- scCNVis::makeSCHeatmap(obj)
```

The output of the plotting functions is a list with two items. The first is the
`Plot` and the second is the `PlotData` that was used to make the plot in case
you wish to create your own version of the plot or inspect the data more thoroughly.
If you wish to simply export the plot, it can be done using the wrapper function:

```
scCNVis::saveCustomPlot(heatmap.result, "~/path/to/heatmap.png")
```


