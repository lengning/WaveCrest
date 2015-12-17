#' @title Identify additional dynamic genes based on recovered order
#' @usage WaveCrestIden(Data, Order, Ndg=5, MeanLOD=10)
#' N=20000, NCThre=1000)
#' @param Data gene-by-sample matrix or isoform-by-sample matrix. It should be rescaled to values bwteen
#' [-1,1].
#' @param Order A numeric vector indicates the reconstructed cell order.
#' @param Ndg degree of polynomial.
#' @param MeanLOD genes whose mean is less than MeanLOD are not considered.
#' @return This function may be used to identify additional dynamic genes based on the recovered order.
#' @examples aa <- sin(seq(0,1,.1))
#' bb <- sin(seq(0.5,1.5,.1))
#' cc <- sin(seq(0.9,1.9,.1))
#' dd <- sin(seq(1.2,2.2,.1))
#' res <- WaveCrestENI(c("aa","bb"), rbind(aa,bb,cc,dd), N=50,
#' Conditions = factor(c(rep(1,5),rep(2,6))))
#' res.more <- WaveCrestIden(rbind(cc,dd), res, MeanLOD=0)
#' @author Ning Leng


WaveCrestIden <- function (Data, Order, Ndg=5, MeanLOD=10){
	  expect_is(Data, "matrix")
		data.use <- Data[which(rowMeans(Data)>MeanLOD),]
    if(length(data.use)==0)stop("no gene passed the mean lower limit of detection threhold")		

		data.sc <- Rescale(data.use[,Order])

    fitall.all <- sapply(1:nrow(data.sc),function(j){
              tt=lm(data.sc[j,]~poly(1:ncol(data.sc),Ndg))
              t2=mean(tt$residuals^2)
                    t2  })
    names(fitall.all) <- rownames(data.sc)
    fitall.all.sort <- sort(fitall.all, decreasing=F)
    res <- fitall.all.sort
}
