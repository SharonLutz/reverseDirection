\name{reverseDirection}
\alias{reverseDirection}
\title{reverseDirection}
\description{
Runs the MR Steiger approach on simulated data in the presence or absence of unmeasured confounding of the SNP-mediator-outcome relationship and pleiotropy.
}
\usage{
reverseDirection(nSim = 1000, n = 100, MAF = 0.5, gamma0 = 0, gammaG = 0.2, varX = 1, measurementError = F, delta0 = 0, deltaX = 1, varME = 1, beta0 = 0, betaX = seq(from = 0, to = 1, length.out = 4), pleiotropy = F, betaG = 1, varY = 0.2, sig.level = 0.05, SEED = 1, plot.pdf = T, plot.name = "reverseDirection.pdf")
}
\arguments{
  \item{nSim}{
the number of simulations
}
  \item{n}{
the sample size
}
  \item{MAF}{
the minor allele frequency of the SNP
}
  \item{gamma0}{
the intercept for the linear association of the SNP G with the exposure X
}
  \item{gammaG}{
the slope for the linear association of the SNP G with the exposure X
}
  \item{varX}{
the variance of the exposure X which is generated from a normal distribution
}
  \item{measurementError}{
if measurementError=T then the exposure is generated with measurement error
}
  \item{delta0}{
the intercept for the linear association of the exposure with measurement error 
}
  \item{deltaX}{
the slope for the linear association of the exposure with measurement error
}
  \item{varME}{
the variance of the exposure X with measurement error which is generated from a normal distribution
}
  \item{beta0}{
the intercept for the linear association of the exposure X with the outcome Y
}
  \item{betaX}{
the slope for the linear association of the exposure X with the outcome Y
}
  \item{pleiotropy}{
if pleiotropy=T then the outcome Y is generated with a direct effect of the SNP G on Y
}
  \item{betaG}{
the slope for the linear association of the SNP G with the outcome Y
}
  \item{varY}{
the variance of the outcome Y which is generated from a normal distribution
}
  \item{sig.level}{
the significance level, default=0.05
}
  \item{SEED}{
the seed for the random number generator
}
  \item{plot.pdf}{
if plot.pdf=T then a plot of the percent of simulations where case 1-3 is accepted. If plot.pdf=F then no plot will be created.
}
  \item{plot.name}{
specifies the name of the plot that is created
}
}
\details{
This function outputs the percent of simulations where case 1 (X->Y), case 2 (X<-Y), or case 3 (inconclusive) is detected by using the MR Steiger approach. The percent of simulations where the Steiger p-value and MR p-value is less than alpha is given (pSteiger and pMR, respectively). The correlation between the SNP G and the exposure X (corGX), correlation between the SNP G and the outcome Y (corGY) and the correlation between the exposure X and the outcome Y (corXY).
}
\references{
Hemani G, Tilling K, Davey Smith G (2017) Orienting the causal relationship between imprecisely measured traits using GWAS summary data. PLOS Genetics 13(11): e1007081.
}
\author{
Sharon Lutz
}
\examples{
reverseDirection(nSim=100)
}
\keyword{reverseDirection }
\keyword{MR}
