# WaveCrest


WaveCrest is a statistical approach to
reconstruct gene expression trajectory in single cell RNA-seq experiments with ordered conditions.
WaveCrest contains two modules - the first module implements an extended nearest insertion (ENI) algorithm and
a 2-opt algorithm that
search for optimal cell orders, and the second module implements a 'fishing' algorithm
that can be used to identify additional dynamic genes following the recovered order.

## 1. Installation
To run the WaveCrest graphical user interface (GUI), it requires the following packages: shiny, shinyFiles, colourpicker, WaveCrest

R version ≥ 3.0.2 is needed. For mac user, make sure whether xcode is installed.

To install the shiny packages, in R run:

> install.packages("shiny")

> install.packages("shinyFiles")

> install.packages("colourpicker")

WaveCrest R package and its vignette could be found at https://github.com/lengning/WaveCrest/tree/master/package

To install WaveCrest, in R run: 

> install.packages("devtools")

> library(devtools)

> install_github("lengning/WaveCrest/package/WaveCrest")

Or install locally.

R users may run WaveCrest following the vignette here https://github.com/lengning/WaveCrest/blob/master/package/WaveCrestVig_v1.pdf and skip the installation of shiny and shinyFiles. 

### Run the app
To launch WaveCrest GUI, in R run:

> library(shiny)

> library(shinyFiles)

> library(colourpicker)

> runGitHub('lengning/WaveCrest')

![Screenshot](https://github.com/lengning/WaveCrest/blob/master/wavecrestscreen.png)

## 2. Input files

The first input file should be the expression matrix, with genes in rows and cells as columns. 
Currently, only takes csv files or tab delimited file are accepted.
The input file will be treated as a tab delimited file if the suffix is not '.csv'.


The second input file is the Condition vector. The conditions could be time points, spatial positions, etc. 
It can be csv or tab delimited file. The file should contain 1 column. Each element of the should represent the corresponding condition that each cell belongs to, it should match exactly the order of columns in the expression matrix and be the same length. If condition input file is missing, all cells are considered to be from one condition.

The third input file is the Marker vector. It can be csv or tab delimited file. The file should contain
1 column of gene names. These genes will be used to do the reordering in WaveCrest. If the marker input file is missing, all genes will considered as markers of interest. If a marker gene is not included in the expression matrix, the marker will be excluded from the analysis.

### Example files
An example input file **exData.csv**, **Condition.csv**, and **Markers.csv** can be found at https://github.com/lengning/WaveCrest/tree/master/example_data   
- Expression matrix contains 200 genes and 120 cells.
- Condition vector indicating there are 30 cells for each time point.
- Marker vector contains a list of 8 marker genes.



## 3. Customize options

By default, the WaveCrest GUI will first recover the cell order based on the markers of interest. In the recovered order, cells from different conditions (time points) are not allowed to be mixed.
If specified, the WaveCrest GUI will further detect additional genes with dynamic profile following the recovered cell order.  

- The number of iteration for 2-opt: Default is 20,000. Increasing the number of iterations may improve the ordering results, but will result in a longer run time.

- Need normalization?: If Yes, normalization using median-by-ratio will be performed prior to the WaveCrest run. If the input matrix is already normalized (e.g. by median-by-ratio normalization or TMM), this option should be disabled by selecting No. In addition, if the input expression matrix only contains a small subset of genes, it is suggested to first perform the normalization using all genes before taking the subset.

- Use log data for input into WaveCrest?: If Yes, data will be logged for the cell reordering. This may alleviate effects due to outliers and produce cleaner results.

- Identify additional dynamic genes based on the recovered order?: If Yes, additional genes that have dynamic profiles similar to the recovered cell order will be identified.

- What type of trend do you expect?: If we assume the target temporal pattern of each marker is monotone increasing (decreasing) from the first cell to the last cell, **Linear** should be used. If the target temporal pattern is expected to follow a quadratic / cubic / quartic polynomial,  **Quadratic / Cubic / Quartic** should be selected. An example quartic form may be bi-modal expression over time. It is important to consider what trend of expression is expected. This will be used to recover the cell order as well as in identifying additional genes.

- Set seed (random number generator): This allows the user to reproduce the results. The same seed number will produce the same results (given the other input is exactly the same).  

- Plot in log scale?: Whether to plot the gene expression on the log scale.


- Plot key markers following recovered cell order?: If Yes, expression plots of input markers following the recovered cell order will be generated. 

- Plot additional dynamic genes following recovered cell order?: If Yes, expression plots of additional dynamic genes following the recovered cell order will be generated. 

- Add trendline to dynamic gene plots?: If Yes, scatterplots will display a fitted polynomial line equal to the degree trend chosen to run WaveCrest (i.e. linear / quadratic / cubic or quartic polynomial).

- Number of additional genes to plot: Number of additional dynamic genes to plot when the previous option is "Yes". If it is not specified, the top 10 genes will be included in the output plots.

- Reverse the recovered cell order: The starting point for WaveCrest is random. Instead of rerunning with another seed, if Yes, this option will just reverse the recovered order. Just make sure to reuse the original seed number!

- Plot heatmap of marker genes following the recovered cell order?: If Yes, then a heatmap will be generate containing the genes from the Marker list with columns ordered according to the recovered cell order.
- Colors for heatmap: If the box is checked then the default geen red colors will be used, otherwise uncheck this box and select 3 new colors to represent low, middle, and high expression in the heatmap.

- Cluster key makers genes in heatmap?: If Yes, then genes will be ordered according to hierarchal clustering, otherwise the original marker list order will be used.

- Output directory, will be set as home directory (~/) if it is empty.
- Output file name for input paramaters and version info.
- Output file name for the normalized expression following original cell order.
- Output file name for the normalized expression following recovered cell order.
- Output file name for the genes sorted by strength of dynamics. This file will not be generated if the second option ("Identify additional dynamic genes") is "No".
- Output file name for the plots (key markers following recovered order).
- Output file name for the plots (additional genes following recovered order). This file will not be generated if the second option ("Identify additional dynamic genes") is "No".
- Output file name for the heatmap (marker genes following recovered order). This file will not be generated if the option ("Plot heatmap of marker genes") is "No".

## 4. Outputs
Two to five files will be generated:
-	normalized.csv: normalized expression matrix with genes in row and cells in column following original cell order.
-	normalized_ENI.csv: normalized expression matrix with genes in row and cells in column following recovered cell order.

- PlotMarkers.pdf: This file will be generated only when the user chooses to plot key markers. In each plot, x-axis shows cells following recovered order and y-axis shows normalized expression. 

-	genes_by_dynamic.csv: This file will be generated only when the user chooses to identify additional dynamic genes. Genes (non-marker) are listed with their corresponding mean squared errors (MSEs). Genes are sorted by MSE.

- PlotDynamic.pdf: This file will be generated only when the user chooses to identify additional dynamic genes and choose to plot additional dynamic genes. In each plot, x-axis shows cells following recovered order and y-axis shows normalized expression. 
 
- WaveCrest_info.txt: This file contains all input parameters and version info.
 
## Note
The 'create new folder' button in the output folder selection pop-up is disfunctional right now



## License
This project is licensed under the terms of the Apache License 2.0

