#' @title Opt-2 algorithm for non-circular Extended Nearest Insertion
#' @usage Opt2TC(Data,N,Seq,condn, Ndg=3, NCThre=1000)
#' @param Data gene-by-sample matrix or isoform-by-sample matrix.It should be r
#' escaled to gene specific z scores.
#' @param condn a numerical vector indicates conditions
#' @param N,NCThre The 2-opt algorithm will stop if N iterations has been performed or if the optimal order 
#' remains unchanged for over NCThre iterations.
#' @param Ndg degree of polynomial.
#' @param Seq a vector indicates the sample order obtained from the ENI.
#' @return This function performs the 
#' the 2-opt algorithm to improve the optimal order searching of the Extended Nearest Insertion (ENI).
#' In each iteration, the function will randomly choose two points (samples), then flip the samples
#' between these two points. The new order will be adapted if it provides smaller PR MSE.
#' The output returns the optimal order and its PR MSE. 
#' Note that the reordering happens within condition (time point). Cells from
#' different conditions won't be mixed unless the cell is placed in the
#' coundary of two conditions.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' res <- ImpTC(rbind(aa,bb,cc), Seq=NULL, condn = c(rep(1,5),rep(2,6)))
#' res2 <- Opt2TC(rbind(aa,bb,cc), N=50, Seq=res, condn= c(rep(1,5),rep(2,6)))
#' @author Ning Leng

Opt2TC <- function(Data,N,Seq,condn,Ndg=3, NCThre=1000){
	message("2-opt: ")
	Ncol <- ncol(Data)
	conduse <- condn[Seq]
	SeqIter <- Seq
	StatIter <- PipeRCDF(Data[,Seq],Ndg)
	nc=0
	condl <- unique(conduse)
	st <- sapply(1:length(condl),function(i)which(conduse==condl[i])[1])
	ed <- Ncol
	if(length(condl)>1)ed <- c(sapply(1:length(condl),function(i)which(conduse==condl[i])[1]-1)[-1],Ncol)
	for(i in 1:N){
		choose.l <- sample(condl,1)
		Choose <-sample(st[choose.l]:ed[choose.l],2)
		Seq0 <- SeqIter
		Seq0[Choose[1]:Choose[2]] <- SeqIter[Choose[2]:Choose[1]]
		Stat <- PipeRCDF(Data[,Seq0],Ndg)
		nc <- nc+1
		if(Stat<StatIter){
			SeqIter=Seq0
		  StatIter=Stat
			  nc=0
			  message("update ",i)
		}
		if(nc > NCThre){
			message("final update:",i)
			break
		}}
Out=list(SeqIter, StatIter)
}

