# WaveCrest


WaveCrest is a statistical approach to
reconstruct gene expression trajectory in single cell RNA-seq experiments with ordered conditions.
WaveCrest contains two modules - the first module implements an extended nearest insertion (ENI) algorithm that
searches for optimal cell orders, and the second module implements a spline fitting module
that can be used to identify additional dynamic genes.

## 1. Installation
The WaveCrest GUI requires the following packages: shiny, shinyFiles, WaveCrest

To install the shiny packages, in R run:

> install.packages("shiny")

> install.packages("shinyFiles")

WaveCrest R package and its vignette could be found at https://github.com/lengning/WaveCrest/tree/master/package

To install WaveCrest, in R run: 

> install.packages("devtools")

> library(devtools)

> install_github("lengning/WaveCrest/package/WaveCrest")

Or install locally.

### Run the app
In R, run:

> library(shiny)

> runGitHub('lengning/WaveCrest')

![Screenshot](https://github.com/lengning/WaveCrest/blob/master/figs/wavecrestscreen.png)

## 2. Input files

The first input file should be the expression matrix. 
Rows are the genes and columns are the cells.
Currently the program only takes csv files or tab delimited file.
The input file will be treated as a tab delimited file if the suffix is not '.csv'.


The second input file is the condition (time points) vector. It could be csv or tab delimited file. The file should contain
1 column. The i th component represents the condition that cell i belongs to. The length of the condition vector should be the same as the number of cells in the first input file. Two or more conditions are expected. If condition input file is missing, all cells are considered to be one condition.

The third input file is the marker vector. It could be csv or tab delimited file. The file should contain
1 column. If marker input file is missing, all genes will considered as markers. If a marker is not included in the expression matrix, the marker will be excluded for the analysis.

### Example files
An example input file **exData.csv**, **Condition.csv**, and **Markers.csv** could be found at https://github.com/lengning/WaveCrest/   
- Expression matrix contains 200 genes and 120 cells 
- Condition vector shows there are 30 cells in each time points
- Marker vecter contains the list of 8 markers

## 3. Customize options
- The number of iteration for ENI: Default is 20000. 
-	Identify additional dynamic genes based on the recovered order?: If Yes, users can 'rescue' gene that were not considered as markers
- What type of trend do you expect?: If we assume the target temporal pattern of each marker is monotone increasing (decreasing) from the first cell to the last cell, **Linear** should be used. If the target temporal pattern is expected to follow a quadratic / cubic / quartic polynomial,  **Quadratic / Cubic / Quartic** should be selected. An example quartic form may be bi-modal expression over time. 
-	Set seed (random number generation): The users can reproduce the results by setting the same seed for different runs.
- Plot key markers following recovered cell order?: If Yes, plot of input markers following recovered cell order will be generated. 
- Plot additional dynamic genes following recovered cell order?: If Yes, plot of additional dynamic genes following recovered cell order will be generated. 
-	Number of additional genes to plot: Number of additional dynamic genes to plot when the previous option is "Yes". If it is not specified, top 10 genes will be included in the output plots.
- Plot in log scale?: Whether take a log scale in plot.
- Output directory, will be set as home directory (~/) if it is empty.
- Output file name for the normalized expression following original cell order.
- Output file name for the normalized expression following recovered cell order.
-	Output file name for the genes sorted by strength of dynamics. This file will not be generate if the second option ("Identify additional dynamic genes") is "No".
-	Output file name for the plots (key markers following recovered order).
-	Output file name for the plots (additional genes following recovered order). This file will not be generate if the second option ("Identify additional dynamic genes") is "No".

## 4. Outputs
Two to five files will be generated:
-	normalized.csv: normalized expression matrix with genes in row and cells in column following original cell order.
-	normalized_ENI.csv: normalized expression matrix with genes in row and cells in column following recovered cell order.

- PlotMarkers.pdf: This file will be generated only when the user chooses to plot key markers. In each plot, x-axis shows cells by recovered order and y-axis shows expression. 

-	genes_by_dynamic.csv: This file will be generated only when the user chooses to identify additional dynamic genes. Genes (non-marker) are listed with their corresponding MSEs. Genes are sorted by MSE.

- PlotDynamic.pdf: This file will be generated only when the user chooses to identify additional dynamic genes and choose to plot additional dynamic genes. In each plot, x-axis shows cells by recovered order and y-axis shows expression. 
 
## Note
The 'create new folder' button in the output folder selection pop-up is disfunctional right now




