\name{WaveCrestExData}
\alias{WaveCrestExData}
\docType{data}
\title{
	Simulated gene level data set with 4 ordered conditions
}
\description{
'WaveCrestExData' gives the gene level simulated data with 4 ordered 
conditions (time points), 30 cells for each condition.
A total of 200 genes were simulated.
Among these 200 genes, 20 genes were simulated to have dynamic profile over the time course data.	
The other 180 genes were simulated as random noise from normal distribution.
The cell order within each condition were perturbed.
}
\format{
	WaveCrestExData is a matrix with 200 genes (rows) and 120 cells (columns). 
}
\seealso{
}
\examples{
data(WaveCrestExData)
str(WaveCrestExData)
}
\keyword{datasets}




