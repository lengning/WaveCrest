library(shiny)
library(shinyFiles)
#library(gdata)
options(shiny.maxRequestSize=500*1024^2) 
# Define UI for slider demo application
shinyUI(pageWithSidebar(
  #  Application title
  headerPanel("Wave-Crest"),
  
  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(width=12,height=20,
               # file
               fileInput("filename", label = "Data file input (support .csv, .txt, .tab)"),
               
               # grouping vector
               fileInput("ConditionVector", label = "Condition vector \n file name (e.g. collection time. support .csv, .txt, .tab)"),
               
               # List of key markers 
               fileInput("Markers", label = "List of key markers \n file name (support .csv, .txt, .tab)"),
               
               column(4,
                      # Num permutation
                      numericInput("Permu",
                                   label = "The number of iteration for 2-opt",
                                   value = 20000),
                      
                      # Normalization
                      radioButtons("Norm_buttons",
                                   label = "Do you need to normalize data?",
                                   choices = list("Yes" = 1,
                                                  "No" = 2),
                                   selected = 1),
				   #Log data for input
                   radioButtons("logDataIn",
                                label = "Use log data for input into WaveCrest?",
                                choices = list("Yes" = 1,
                                               "No" = 2),
                                selected = 1),
								
                      # Identify additional genes
                      radioButtons("Iden_buttons",
                                   label = "Identify additional dynamic genes based on the recovered order ('fishing')?",
                                   choices = list("Yes" = 1,
                                                  "No" = 2),
                                   selected = 1),
                      
                      # Identify additional genes
                      radioButtons("DF_buttons",
                                   label = "What type of trend do you expect?",
                                   choices = list("Linear" = 1,
                                                  "Quadratic" = 2,
                                                  "Cubic"=3,
                                                  "Quartic"=4),
                                   selected = 3),
                      # set seed
                      numericInput("Seed", 
                                   label = "Set seed (for random number generator)", 
                                   value = 1),
                      # plot log-exp or not
                      radioButtons("log_whether",
                                   label = "Plot in log scale?",
                                   choices = list("No" = 1,
                                                  "log2(expression + 1)" = 2),
                                   selected = 1),
                      
                      radioButtons("MarkerPlot_buttons",
                                   label = "Plot key markers following recovered cell order?",
                                   choices = list("Yes" = 1,
                                                  "No" = 2),
                                   selected = 1)
                      
               ),
               
               column(width=4,
                     
                      
                      radioButtons("AddPlot_buttons",
                                   label = "Plot additional dynamic genes following recovered cell order?",
                                   choices = list("Yes" = 1,
                                                  "No" = 2),
                                   selected = 1),        				
                      # num genes to plot
                      textInput("PlotNum", 
                                label = "Number of additional genes to plot (if not specified, top 10 genes will be plotted)", 
                                value = ""),
                      
                   
								   
								   
                      radioButtons("AddHeatMap_buttons",
                                    label = "Plot heatmap of key marker genes following recovered cell order?",
                                    choices = list("Yes" = 1,
                                                   "No" = 2),
                                    selected = 1),  
		
									
                
                      tags$div(tags$b("Select colors for the heatmap :")),
                      checkboxInput("UseCols", "Use default green to red heatmap ? ", TRUE),
                      fluidPage(colourInput("col1", "Low values: ", "green")),
                      fluidPage(colourInput("col2", "Middle values:", "black")),
                      fluidPage(colourInput("col3", "High values:", "red")),
                      
                      radioButtons("clusterRow",
                                   label = "Cluster key marker genes in heatmap (otherwise use original order) ? ",
                                   choices = list("Yes" = 1,
                                                  "No" = 2),
                                   selected = 1)
                      
                
               ),
               
               column(4,
                      # output dir
                      tags$div(tags$b("Please select a folder for output :")),
                      
                      shinyDirButton('Outdir', label ='Select Output Folder', title = 'Please select a folder'),
                      tags$br(),
                      tags$br(),
                      
                      # plot name
                      textInput("InfoFileName", 
                                label = "Export file name - input parameters and version info", 
                                value = "WaveCrest_info"),
               
                      # export normalzied matrix (original order)
                      textInput("exNormFileName", 
                                label = "Export file name - normalized expression matrix (following original cell order)", 
                                value = "normalized"),
                      
                      # export normalzied matrix (recovered order)
                      textInput("exENINormFileName", 
                                label = "Export file name - normalized expression matrix (following recovered cell order)", 
                                value = "normalized_ENI"),                
                      
                      # export gene list with p-value (IdenRes)
                      textInput("exGListFileName", 
                                label = "Export file name - genes sorted by fishing MSE", 
                                value = "genes_by_dynamic"),
                      # plot name
                      textInput("exMarkerPlotFileName", 
                                label = "Export file name for the plots? (key markers following recovered order)", 
                                value = "PlotMarkers"),
                      # plot name
                      textInput("exDynamicPlotFileName", 
                                label = "Export file name for the plots? (additional genes following recovered order)", 
                                value = "PlotDynamic"),
								
                      # plot name
                      textInput("exHeatMapPlotFileName", 
                                label = "Export file name for the plots? (heatmap of key markers in recovered order)", 
                                value = "PlotHeatMap")
         
                      
                      ),
               br(),
               br(),
               actionButton("Submit","Submit for processing")
  ),
  
  # Show a table summarizing the values entered
  mainPanel(
    h4(textOutput("print0")),
    #tableOutput("values")
    dataTableOutput("tab")
  )
))
