\name{reverseDirection}
\alias{reverseDirection}
\title{reverseDirection}
\description{
Runs the MR Steiger approach on simulated data in the presence or absence of unmeasured confounding of the SNP-mediator-outcome relationship.
}
\usage{
reverseDirection(nSim = 1000, n = 100, nSNP = 1, MAF = rep(0.5, 1), gamma0 = 0, gammaX = rep(0.2, 1), varM = 1, beta0 = 0, betaM = seq(from = 0, to = 1, length.out = 4), varY = 0.2, delta0 = 0, deltaX = rep(1, 1), varU = 1, gammaU = 0, betaU = 1, Uconfounder = T, alpha = 0.05, SEED = 1, plot.pdf = T, plot.name = "reverseDirection.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nSim}{
the number of simulations
}
  \item{n}{
the sample size
}
  \item{nSNP}{
the number of SNPs associated with the mediator M
}
  \item{MAF}{
a vector of the minor allele frequency of the SNPs
}
  \item{gamma0}{
the intercept for the linear association of the SNPs with the mediator M
}
  \item{gammaX}{
the slope for the linear association of the SNPs with the mediator M
}
  \item{varM}{
the variance of the mediator M which is generated from a normal distribution
}
  \item{beta0}{
the intercept for the linear association of the mediator M with the outcome Y
}
  \item{betaM}{
the slope for the linear association of the mediator M with the outcome Y
}
  \item{varY}{
the variance of the outcome Y which is generated from a normal distribution
}
  \item{delta0}{
the intercept for the linear association of the SNPs with the unmeasured confounder U
}
  \item{deltaX}{
the slope for the linear association of the SNPs with the unmeasured confounder U
}
  \item{varU}{
the variance of the unmeasured confouncer U which is generated from a normal distribution
}
  \item{gammaU}{
the slope for the linear association of the SNPs with the unmeasured confounder U
}
  \item{betaU}{
the slope for the linear association of the unmeasured confounder U with the outcome Y
}
  \item{Uconfounder}{
if Uconfounder=T then an unmeasured confounder is genreated if Uconfouder=F then no unmeasured confounder is generated. The unmeausred confounder is generated from a normal distribution based on the parameters defined above.
}
  \item{alpha}{
the significance level
}
  \item{SEED}{
the seed for the random number generator
}
  \item{plot.pdf}{
if plot.pdf=T then a plot of the percent of simulations where the correct direction is detected will be generated. If plot.pdf=F then no plot will be created.
}
  \item{plot.name}{
specifies the name of the plot that is created
}
}
\value{
This function outputs the percent of simualtions where the correct direction is detected  between the mediator M and outcome Y using the MR Steiger approach either not adjusting for measurment error or adjusting for measurment error (CorrectDirection and CorrectDirectionAdj, respectively). The percent of simulations where the Steiger p-value is less than alpha without or with adjusting for measurment error (SteigerTest and SteigerTestAdj, respectively) The average sensitivity ratio is given by SensitivityRatio. The correlation between the first SNP and the mediator M (corX1M), correlation between the first SNP and the outcome Y (corX1Y) and the correlation between the mediator M and the outcome Y (corMY).
}
\references{
Hemani G, Tilling K, Davey Smith G (2017) Orienting the causal relationship between imprecisely measured traits using GWAS summary data. PLOS Genetics 13(11): e1007081.
}
\author{
Sharon Lutz
}
\examples{
reverseDirection(nSim=1000)
}
\keyword{reverseDirection }
\keyword{MR}
