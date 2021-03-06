#' @title Search for the optimal sample order for experiments with ordered conditions
#' @usage WaveCrestENI(GeneList, Data, Conditions, Ndg=3,
#' N=20000, NCThre=1000)
#' @param Data gene-by-sample matrix or isoform-by-sample matrix. It should be rescaled to values bwteen
#' [-1,1].
#' @param GeneList A vector that contains genes to use for reordering.
#' @param Conditions A factor indicates the condition (time/spatial point) which
#' each sample belongs to. The levels of the factor should be sorted
#' by the time-course / spatial order.
#' @param Ndg degree of polynomial.
#' @param N,NCThre The 2-opt algorithm will stop if N iterations has been performed or if the optimal order 
#' remains unchanged for over NCThre iterations.
#' @return This function performs the extended nearest insertion (ENI) and 2-opt algorithm.
#' The ENI algorithm will be applied to
#' search for the optimal sample order which minimizes the MSE of 
#' polynomial regression (PR). 
#' This function will call PipeRCDF() function, which fits 
#' PR to expression of each gene/isoform within the gene list. 
#' The aggregated MSE of a fit is defined as the 
#' summation of the MSEs of all genes/isoforms considered here.
#' The output of PipeRCDF() returns the optimal order which provides the smallest PR MSE.
#' The 2-opt algorithm will then be applied to improve the optimal order searching of the ENI.
#' In each iteration, the 2-opt algorithm will randomly choose two points (samples), then flip the samples
#' between these two points. The new order will be adapted if it provides smaller PR MSE.
#' The output returns the optimal order. 
#' Note that the reordering happens within condition (time point). Cells from
#' different conditions won't be mixed unless the cell is placed in the
#' coundary of two conditions.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' dd <- sin(seq(1.2,2.2,.1))
#' res <- WaveCrestENI(c("aa","bb"), rbind(aa,bb,cc,dd), N=50,
#' Conditions = factor(c(rep(1,5),rep(2,6))))
#' @author Ning Leng



WaveCrestENI <- function (GeneList, Data, Conditions, Ndg=3, N=20000, NCThre=1000, Seed=1000){
	expect_is(Data, "matrix")
	set.seed(Seed)
	start.order <- sample(1:ncol(Data),ncol(Data))
	if(!is.factor(Conditions))Conditions <- factor(Conditions, levels = unique(Conditions))
	cond.num <- as.numeric(Conditions)
	
	if(length(setdiff(GeneList,rownames(Data)))>0)stop("some genes in GeneList are not in Data!")
	data.use <- Data[GeneList,]
	data.sc <- Rescale(data.use)
	res.imp <- ImpTC(data.sc, start.order,cond.num, Ndg=Ndg)
	res.out <- Opt2TC(data.sc, N, res.imp, cond.num, Ndg=Ndg, NCThre=NCThre)
	res.out[[1]]
}
