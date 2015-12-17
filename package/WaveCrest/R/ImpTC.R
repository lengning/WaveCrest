#' @title Non-circular Extended Nearest Insertion
#' @usage ImpTC(Data, Seq, condn, Ndg=3)
#' @param Data gene-by-sample matrix or isoform-by-sample matrix.It should be rescaled to gene specific z scores.
#' @param Ndg degree of polynomial.
#' @param Seq NULL or a vector indicates the sample order.
#' if specified, the samples will be first reordered by this vector.
#' @param condn a numerical vector indicates conditions
#' @return This function performs the extended nearest insertion (ENI).
#' The ENI algorithm searchs for the optimal sample order which minimizes the MSE of 
#' the polynomial regression  (PR). 
#' This function will call PipeRCDF() function, which fits 
#' PR to each row of the data. 
#' The aggregated MSE of a fit is defined as the 
#' summation of the MSEs of all genes/isoforms considered here.
#' The output returns the optimal order which provides the smallest PR MSE.
#' Note that the reordering happens within condition (time point). Cells from
#' different conditions won't be mixed unless the cell is placed in the
#' coundary of two conditions.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' res <- ImpTC(rbind(aa,bb,cc), Seq=NULL, condn = c(rep(1,5),rep(2,6)))
#' @author Ning Leng
ImpTC <- function(Data, Seq, condn, Ndg=3){
	if(is.null(Seq))Seq=1:ncol(Data)
  expect_is(Seq, "integer")
	conduse <- condn[Seq]
	
	Od <- Seq[1:3]
	Od <- Od[order(conduse[1:3])] # sort by conditions
	condnow <- condn[Od]
	Ncol <- ncol(Data)
	for(i in 4:ncol(Data)){
	tp <- Seq[i]
	tpc <- conduse[i]
	w1 <- which(condnow<tpc)
	w2 <- which(condnow>tpc)
	condleft <- ifelse(length(w1)==0, NA, max(w1))
	condright <- ifelse(length(w2)==0, NA, min(w2))
	#matrix with i columns, i rows
 		
	if(!is.na(condleft) & !is.na(condright))
	mat=t(sapply(condleft:(condright-1),function(j)c(Od[1:j],tp,Od[(j+1):(i-1)])))
	if(is.na(condleft))mat <- matrix(c(tp, Od[1:(i-1)]),nrow=1)
	if(is.na(condright))mat <- matrix(c( Od[1:(i-1)], tp),nrow=1)
	if(is.na(condright) & is.na(condleft)) mat <- rbind(c(tp, Od[1:(i-1)]), c( Od[1:(i-1)], tp))
	#Ndg=ceiling(i/Seg)
	Stat=sapply(1:nrow(mat),function(j)PipeRCDF(Data[,mat[j,]], Ndg))
	Min=which.min(Stat)
	Od=mat[Min,]
	condnow=condn[Od]
	message("insert ", i , "th cell")
	#print(tpc)
	#print(condnow)
	#browser()
	}
Od
}


