# reverseDirection
Examines the MR Steiger approach to detect the directionality between the mediator and outcome through simulation studies.

## Installation
```
# R v3.4 or higher is needed
install.packages("devtools")  # devtools must be installed first
install.packages("MendelianRandomization") 
install.packages("psych")
devtools::install_github("MRCIEU/TwoSampleMR") 

devtools::install_github("SharonLutz/reverseDirection")
```

## Input
First, the number of SNPs (input nSNP) are generated from a binomial distribution for n subjects (input n) for a given minor allele frequency (input vector MAF).

For the SNP Xi for i=1,...,K for K SNPs, the mediator/ exposure M is generated from a normal distribution with the variance (input varM) and the mean as follows:

E\[M \] = &gamma;<sub>o</sub> + &sum; &gamma;<sub>X</sub>  X<sub>i</sub> 

All of these values are inputted by the user (i.e. the intercept gamma0, and the genetic effect size as a vector gammaX).

The outcome Y is generated from a normal distribution with the variance (input varY) and the mean as follows:

E\[Y \] = &beta;<sub>o</sub> +  &beta;<sub>M</sub> M

All of these values are inputted by the user (i.e. the intercept beta0 and the effect of the mediator directly on the outcome as betaM).

If there is pleiotropy (input Uconfounder =T), then an additional covariate U is gernerated from a normal distribution with a variance (input varU) and a mean such that

E\[U \] = &delta;<sub>o</sub> + &sum; &delta;<sub>X</sub>  X<sub>i</sub> 

Then, the mediator is generated as defined above, but the outcome Y is generated from a normal distribution such that
E\[Y \] = &beta;<sub>o</sub> +  &beta;<sub>M</sub> M  +  &beta;<sub>U</sub> U

After the SNPs X, mediator M, and outcome Y are generated, then the reverseDirection function runs the MR Steiger approach to determine if the mediator M causes the outcome Y.

## Output
This function outputs the percent of simualtions where the correct direction is detected between the mediator M and outcome Y using the MR Steiger approach either not adjusting for measurment error or adjusting for measurment error (CorrectDirection and CorrectDirectionAdj, respectively). The percent of simulations where the Steiger p-value is less than alpha without or with adjusting for measurment error (SteigerTest and SteigerTestAdj, respectively) The average sensitivity ratio is given by SensitivityRatio. The correlation between the first SNP and the mediator M (corX1M), correlation between the first SNP and the outcome Y (corX1Y) and the correlation between the mediator M and the outcome Y (corMY).

## Example:
Consider an example with 100 subjects (input n=100) for one SNP (input nSNP = 1) with a MAF of 50% (input MAF=0.5). Consider a pleiotropic effect (input Uconfounder =T). Then, let the mediator M be generated from a normal distribution with a variance of 1 (input varM = 1) and mean such that 
E\[M \] = 0 + 0.4 X
(input gamma0 = 0, gammaX 0.4). The covariate U is generated from a normal distribution with a variance of 1 (input varU = 1) and mean such that 
E\[U \] = 0 + 0.25 X
(input delta0 = 0, deltaX = 0.25). The outcome Y is generated from a normal distribution with a variance of 1 (input varU = 1) and mean such that 
E\[Y \] = 0 + &beta;<sub>M</sub> X +0.25 U
(input beta0 = 0, betaU = 0.25) and &beta;<sub>M</sub> varies from 0 to 1 by 0.25 (input betaM =seq(from = 0, to = 1, by=0.25)). The R code to run this example is given below.

```
library(reverseDirection)

rr<-reverseDirection(nSim = 1000, n = 100, nSNP = 1, MAF = 0.5, gamma0 = 0, gammaX = 0.4, varM = 1, 
beta0 = 0, betaM =seq(from = 0, to = 1, by=0.25) , varY = 0.2, delta0 = 0, deltaX = 0.25, 
varU = 1, gammaU = 0, betaU = 0.25, Uconfounder =T, alpha = 0.05, SEED = 1, 
plot.pdf = T, plot.name = "plotMRdirection.pdf")

round(rr$matrix,2)
```

```
     McausesY McausesYadj CorrectDirection CorrectDirectionAdj SteigerTest SteigerTestAdj
[1,]     0.18        0.18             0.90                0.90        0.18           0.18
[2,]     0.04        0.04             0.78                0.78        0.04           0.04
[3,]     0.00        0.00             0.61                0.61        0.00           0.00
[4,]     0.00        0.00             0.47                0.47        0.00           0.00
[5,]     0.00        0.00             0.44                0.44        0.00           0.00
     SensitivityRatio corX1M corX1Y corMY
[1,]            55.66   0.27   0.09  0.02
[2,]            76.08   0.27   0.20  0.47
[3,]             2.98   0.27   0.25  0.72
[4,]             2.27   0.27   0.27  0.84
[5,]             1.76   0.27   0.28  0.90
```

<img src="plotMRdirection.png" width="400">


