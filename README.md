# WaveCrest


WaveCrest is a statistical approach to
reconstruct gene expression trajectory in single cell RNA-seq experiments with ordered conditions.
WaveCrest contains two modules - the first module implements an extended nearest insertion (ENI) algorithm and
a 2-opt algorithm that
search for optimal cell orders, and the second module implements a 'fishing' algorithm
that can be used to identify additional dynamic genes following the recovered order.

## 1. Installation
To run the WaveCrest graphical user interface (GUI), it requires the following packages: shiny, shinyFiles, WaveCrest

To install the shiny packages, in R run:

> install.packages("shiny")

> install.packages("shinyFiles")

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

> runGitHub('lengning/WaveCrest')

![Screenshot](https://github.com/lengning/WaveCrest/blob/master/figs/wavecrestscreen.png)

## 2. Input files

The first input file should be the expression matrix. 
Rows are the genes and columns are the cells.
Currently the program only takes csv files or tab delimited file.
The input file will be treated as a tab delimited file if the suffix is not '.csv'.


The second input file is the condition vector. The conditions could be time points, spatial positions, etc. 
It could be csv or tab delimited file. The file should contain
1 column. The i th component represents the condition that cell i belongs to. The length of the condition vector should be the same as the number of cells in the first input file. Two or more conditions are expected. If condition input file is missing, all cells are considered to be from one condition.

The third input file is the marker vector. It could be csv or tab delimited file. The file should contain
1 column, elements are the gene names.
If marker input file is missing, all genes will considered as markers of interest. If a marker is not included in the expression matrix, the marker will be excluded for the analysis.

### Example files
An example input file **exData.csv**, **Condition.csv**, and **Markers.csv** could be found at https://github.com/lengning/WaveCrest/tree/master/example_data   
- Expression matrix contains 200 genes and 120 cells 
- Condition vector shows there are 30 cells in each time points
- Marker vecter contains the list of 8 markers



## 3. Customize options

By default, the WaveCrest GUI will first recover the cell order based on the markers of interest. In the recovered order, cells from different conditions (time points) are not allowed to be mixed.
If specified, the WaveCrest GUI will further detect additional genes with dynamic profile following the recovered cell order.  

- Need normalization? If Yes, normalization will be performed prior to WaveCrest run. If the input matrix is normalized (e.g. by median-by-ratio normalization or TMM), this option should be disabled. In addition, if the input expression matrix only contains a small set of genes, it is suggested to normalize using all genes first before taking the subset.
- The number of iteration for 2-opt: Default is 20000. Increasing the number of iteration may improve the ordering results, but will result in longer run time.
-	Identify additional dynamic genes based on the recovered order?: If Yes, users can identify additional genes that have dynamic profile following the recovered cell order.
- What type of trend do you expect?: If we assume the target temporal pattern of each marker is monotone increasing (decreasing) from the first cell to the last cell, **Linear** should be used. If the target temporal pattern is expected to follow a quadratic / cubic / quartic polynomial,  **Quadratic / Cubic / Quartic** should be selected. An example quartic form may be bi-modal expression over time. 
-	Set seed (random number generator): The users can reproduce the results by setting the same seed for different runs.
- Plot key markers following recovered cell order?: If Yes, expression plots of input markers following recovered cell order will be generated. 
- Plot additional dynamic genes following recovered cell order?: If Yes, expression plots of additional dynamic genes following recovered cell order will be generated. 
-	Number of additional genes to plot: Number of additional dynamic genes to plot when the previous option is "Yes". If it is not specified, top 10 genes will be included in the output plots.
- Plot in log scale?: Whether plot the expressions in log scale.
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

- PlotMarkers.pdf: This file will be generated only when the user chooses to plot key markers. In each plot, x-axis shows cells following recovered order and y-axis shows normalized expression. 

-	genes_by_dynamic.csv: This file will be generated only when the user chooses to identify additional dynamic genes. Genes (non-marker) are listed with their corresponding mean squared errors (MSEs). Genes are sorted by MSE.

- PlotDynamic.pdf: This file will be generated only when the user chooses to identify additional dynamic genes and choose to plot additional dynamic genes. In each plot, x-axis shows cells following recovered order and y-axis shows normalized expression. 
 
## Note
The 'create new folder' button in the output folder selection pop-up is disfunctional right now




