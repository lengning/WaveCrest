#' @title Summarize residuals of polynomial fittings among multiple genes for non-circular Extended Nearest Insertion
#' @usage PipeRCDF(Data, Ndg=20)
#' @param Data gene-by-sample matrix
#' @param Ndg degree in polynomial fitting
#' @return The function will fit polynomial regression  (SPR)
#' to each row of the data. 
#' The aggregated MSE of a fit is defined as the 
#' summation of the MSEs of all genes/isoforms considered here.
#' The output returns the aggregated MSE. 
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' res <- PipeRCDF(rbind(aa,bb,cc), Ndg=3)

#' @author Ning Leng

PipeRCDF <- function(Data,Ndg=20){
	t1 <- PipeR(Data, Ndg)
	EC <- quantile(ecdf(t1), 1:100/100)
	out <- sum(EC)
}

