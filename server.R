library(shiny)
library(shinyFiles)
library(WaveCrest)
#library(EBSeq)
library(colourpicker)

# Define server logic for slider examples
shinyServer(function(input, output, session) {
  volumes <- c('home'="~")
  shinyDirChoose(input, 'Outdir', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$Dir <- renderPrint({parseDirPath(volumes, input$Outdir)})
  
  wcinfo = sessionInfo(package="WaveCrest")
  wcinfo_print = wcinfo$otherPkgs
  
  In <- reactive({
    print(input$Outdir)
    outdir <- paste0("~", do.call("file.path", input$Outdir[[1]]), "/")
    print(outdir)
    
    the.file <- input$filename$name
    if(is.null(the.file))stop("Please upload data")
    Sep=strsplit(the.file,split="\\.")[[1]]
    if(Sep[length(Sep)]=="csv")a1=read.csv(input$filename$datapath,stringsAsFactors=F,header=TRUE, row.names=1)
    if(Sep[length(Sep)]!="csv")a1=read.table(input$filename$datapath,stringsAsFactors=F,header=TRUE, row.names=1)
    Data=data.matrix(a1)
    
    Group.file <- input$ConditionVector$name
    if(is.null(Group.file))GroupVIn = list(c1=rep(1,ncol(Data)))
    if(!is.null(Group.file)){
      Group.Sep=strsplit(Group.file,split="\\.")[[1]]
      if(Group.Sep[length(Group.Sep)]=="csv")GroupVIn=read.csv(input$ConditionVector$datapath,stringsAsFactors=F,header=F)
      if(Group.Sep[length(Group.Sep)]!="csv")GroupVIn=read.table(input$ConditionVector$datapath,stringsAsFactors=F,header=F, sep="\t")
    }
    GroupV=GroupVIn[[1]]
    if(length(GroupV)!=ncol(Data)) stop("length of the condition vector is not the same as number of cells!")
  
    Marker.file <- input$Markers$name
    if(is.null(Marker.file))MarkerVIn = list(c1=rownames(Data))
    if(!is.null(Marker.file)){
      Marker.Sep=strsplit(Marker.file,split="\\.")[[1]]
      if(Marker.Sep[length(Marker.Sep)]=="csv")MarkerVIn=read.csv(input$Markers$datapath,stringsAsFactors=F,header=F)
      if(Marker.Sep[length(Marker.Sep)]!="csv")MarkerVIn=read.table(input$Markers$datapath,stringsAsFactors=F,header=F, sep="\t")  
    }
    MarkerV=MarkerVIn[[1]]
    
	print(input)
    # Compose data frame
    #input$filename$name
    List <- list(
      Input=the.file,
      GroupFile=Group.file,
      MarkerFile=Marker.file,
      Permu=input$Permu, 
      NormTF = ifelse(input$Norm_buttons=="1",TRUE,FALSE), 
	  LogInTF = ifelse(input$Log_InData=="1",TRUE,FALSE), 
      Cond=factor(GroupV, levels=unique(GroupV)),# follow the order they appeared
      Marker=MarkerV,# follow the order they appeared
      test=ifelse(input$Iden_buttons=="1",TRUE,FALSE), 
      testDF=input$DF_buttons, 
      Seed=input$Seed,
      Dir=outdir, 
      exExpF = paste0(outdir,input$exNormFileName,".csv"),
      exENIExpF = paste0(outdir,input$exENINormFileName,".csv"),		
      exDGF = paste0(outdir,input$exGListFileName,".csv"),
      UseCols = ifelse(input$UseCols=="1",TRUE,FALSE), 
      LowCol = input$col1,
      MidCol = input$col2,
      HighCol = input$col3,
      clusterRowT = ifelse(input$clusterRow=="1",TRUE,FALSE), 
      PlotMarkerTF = ifelse(input$MarkerPlot_buttons=="1",TRUE,FALSE), 
      PlotAddTF = ifelse(input$AddPlot_buttons=="1",TRUE,FALSE), 
	  PlotHeatTF = ifelse(input$AddHeatMap_buttons=="1",TRUE,FALSE),
      whetherLog = ifelse(input$log_whether=="1",FALSE,TRUE),
      PlotN = input$PlotNum,
      MarkerPlotF = paste0(outdir,input$exMarkerPlotFileName,".pdf"),
      DynamicPlotF = paste0(outdir,input$exDynamicPlotFileName,".pdf"),  
	  HeatMapF = paste0(outdir,input$exHeatMapPlotFileName,".pdf"),
      Info = paste0(outdir,input$InfoFileName,".txt")
    )
    if(is.null(Marker.file) & List$test==TRUE) print("Warning: All genes are used as markers")
    if(is.null(Marker.file)) List$test=FALSE
    numdegree = switch(List$testDF,"1"=1,"2"=2,"3"=3,"4"=4)
    if(!List$test)List$PlotAddTF = FALSE
    # normalization     
    if(List$NormTF){
    Sizes <- MedianNorm(Data)
    if(is.na(Sizes[1])){
      Sizes <- MedianNorm(Data, alternative=TRUE)
      message("alternative normalization method is applied")
    }
    Data <- GetNormalizedMat(Data,Sizes)
    }
    
    if(!List$NormTF){
      Data <- Data
    }
	
	DataUse <- Data
	
    # main function
    if(length(which(!List$Marker %in% rownames(Data)))>0) {
      print("Warning: not all provided markers are in data matrix")
      List$Marker = intersect(rownames(Data),List$Marker)
    }
	
	if(List$LogInTF == TRUE) {
		DataUse[which(DataUse < 1)] <- 1 # truncate values less than one
	    DataUse <- log2(DataUse) # log2 of the truncated data
	}
	
	#select random starting order
	set.seed(List$Seed)
    rd <- sample(1:ncol(DataUse),ncol(DataUse)) 
    DataUse.rand <- DataUse[,rd]
	Cond.rand <- List$Cond[rd]
    
	#Get cell orders
    ENIRes <- WaveCrestENI(List$Marker, DataUse.rand, Cond.rand, Ndg = numdegree, N = List$Permu, Seed = List$Seed)	
    print("WaveCrestENI...")
    
	
	ENIRes.Order <- colnames(DataUse.rand)[ENIRes]
	 
	
	#Test for additional genes
    if(List$test){
    DataRemain <- DataUse[setdiff(rownames(DataUse),List$Marker),]
    IdenRes <- WaveCrestIden(DataRemain,ENIRes.Order, Ndg = numdegree)
    IdenRes <- sort(IdenRes, decreasing = TRUE)
    DGlist <- cbind(names(IdenRes),IdenRes)
    colnames(DGlist) <- c("gene", "MSE")
    write.csv(DGlist,file=List$exDGF)
    }
    if(!List$test){
      DGlist=c("nope")
    }
      print("Identify additional genes...")
    write.csv(Data, file=List$exExpF) #write input 
    write.csv(Data[,ENIRes.Order], file=List$exENIExpF) #write input with order
    
    
    if(List$PlotMarkerTF){
      PN <- length(List$Marker)
        pdf(List$MarkerPlotF, height=10,width=15)
        par(mfrow=c(3,3))
        for(i in 1:PN){
          if(List$whetherLog==TRUE) plot(log2(Data[List$Marker[i],ENIRes.Order]+1), 
                                        col=as.numeric(List$Cond), ylab="log2(expression+1)",
                                        main=List$Marker[i])
          if(List$whetherLog==FALSE) plot(Data[List$Marker[i],ENIRes.Order], 
                                        col=as.numeric(List$Cond), ylab="expression",
                                        main=List$Marker[i] )  
        }
        dev.off()
		
    }
    
    
	if(List$PlotHeatTF) {
		if(List$UseCols == TRUE) {
			getcols = greenred(100)
			HeatCols = getcols[c(1:30, seq(31,70,3), 71:100)]
		} else if(List$UseCols == FALSE) {
			    getcols = colorpanel(100, List$LowCol, List$MidCol, List$HighCol)
		    	HeatCols = getcols[c(1:30, seq(31,70,3), 71:100)]
			
		}
    
	  data.reorder <- Data[List$Marker,ENIRes.Order]
	  data.reorder[data.reorder < 1] <- 1
	  data.use.log <- log2(data.reorder)
	  data.use.rs <- Rescale(data.use.log)
	  
	  if(length(List$Marker) < ncol(Data)) {hght = 6; wght = 10};
	  if(length(List$Marker) >= ncol(Data)) {hght = 12; wght = 6}; 	 

	  pdf(List$HeatMapF, height=hght,width=wght)
		
		if(List$clusterRowT == TRUE) {
			heatmap.2(data.use.rs, Colv=FALSE, trace='none', reorderfun=function(d,w)rev(d),col = HeatCols, mar = c(5,8))
        } else {	
			heatmap.2(data.use.rs, Colv=FALSE, trace='none', dendrogram='none', Rowv=FALSE, col = HeatCols, mar = c(5,8))
		  }
     dev.off()
		}

    
    if(List$PlotAddTF){      
      PN <- NULL
      if(List$PlotN!="")PN <- as.numeric(List$PlotN)
      else PN <- 10  
      if(length(PN)>0){
        pdf(List$DynamicPlotF, height=10,width=15)
        par(mfrow=c(3,3))
        for(i in 1:PN){
            if(List$whetherLog==TRUE) plot(log2(Data[names(IdenRes)[i],ENIRes.Order]+1), 
                                            col=as.numeric(List$Cond), ylab="log2(expression+1)",
                                            main=names(IdenRes)[i])
            if(List$whetherLog==FALSE) plot(Data[names(IdenRes)[i],ENIRes.Order], 
                                           col=as.numeric(List$Cond), ylab="expression",
                                           main=names(IdenRes)[i])  
        }
        dev.off()
      }}
    
    ## Sessioninfo & input parameters
    sink(List$Info)
    print(paste0("Package version: ", "WaveCrest_",wcinfo_print$WaveCrest$Version))
    print("Input parameters")
    print(paste0("the number of iteration for 2-opt? ", List$Permu))
    print(paste0("whether normalize data? ", List$NormTF))
    print(paste0("whether identify additional dynamic genes based on the recovered order? ", List$test))
    print(paste0("what type of trend the user expect? ", List$testDF))
    print(paste0("seed?: ",List$Seed))
    print(paste0("whether plot key markers following recovered cell order? ", List$PlotMarkerTF))
    print(paste0("whether plot additional genes following recovered cell order? ", List$PlotAddTF))
    print(paste0("how many additional genes to plot? ", PN))
    print(paste0("whether plot in log scale? ", List$whetherLog))
    print(paste0("colors used in the heatmap: Low Color = ", List$LowCol,", Middle Color = ", List$MidCol, ", High Color = ", List$HighCol))
    print(paste0("whether to cluster key markers in heatmap? ", List$clusterRow))
    sink()
    #sink(file="/tmp/none");sink("/dev/null")
    
    
    List=c(List, list(Sig=DGlist))	
  })   
  
  Act <- eventReactive(input$Submit,{
    In()})
  # Show the values using an HTML table
  output$print0 <- renderText({
    tmp <- Act()
    str(tmp)
    paste("output directory:", tmp$Dir)
  })
  
  output$tab <- renderDataTable({
    tmp <- Act()$Sig
    t1 <- tmp
    print("done")
    t1
  },options = list(lengthManu = c(4,4), pageLength = 20))
  

})