#' @title Rescale the gene/isoform expression matrix 
#' @usage Rescale(Data)
#' @param Data input gene-by-sample matrix or isoform-by-sample matrix
#' @note The output will be a gene-by-sample or isoform-by-sample matrix.
#' For each gene/isoform, the expressions will be scaled to mean 0 var 1
#' @examples Rescale(matrix(rnorm(10), nrow=2))
#' @author Ning Leng


Rescale <- function(Data){
InEC=Data
ECNo0=InEC[which(rowMeans(InEC)>0),]
ECSC=t(apply(ECNo0,1,scale))
rownames(ECSC)=rownames(ECNo0)
colnames(ECSC)=colnames(ECNo0)
ECSC
}
