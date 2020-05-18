# reverseDirection
Examines the MR Steiger approach to detect the directionality between the exposure X and outcome Y through simulation studies.

## Installation
```
install.packages("devtools")  # devtools must be installed first, R v3.4 or higher is needed
install.packages("psych")

devtools::install_github("SharonLutz/reverseDirection")
```

## Input
First, the SNP is generated from a binomial distribution for n subjects (input n) for a given minor allele frequency (input MAF).

For the SNP G, the true exposure (X<sub>true</sub>) is generated from a normal distribution with the variance (input varX) and the mean as follows:

E\[Xtrue \] = &gamma;<sub>o</sub> + &gamma;<sub>G</sub> G

All of these values are inputted by the user (i.e. the intercept gamma0, and the genetic effect size gammaG). If there is no measurement error (input measurementError=F), then X=X<sub>true</sub>. If there is measurement error (input measurementError=T), then the measured exposure X is generated from the true exposure X<sub>true</sub> such that

E\[X \] = &delta;<sub>o</sub> + &delta;<sub>X</sub> X<sub>true</sub>

where &delta;<sub>o</sub> and &delta;<sub>X</sub> are inputted by the user (input delta0, deltaX). The outcome Y is generated from a normal distribution with the variance (input varY) and the mean as follows:

E\[Y \] = &beta;<sub>o</sub> +  &beta;<sub>X</sub> X<sub>true</sub>

if there is no pleiotropy (input pleiotropy=F). If there is pleiotropy (input pleiotropy=T), then the outcome Y is generated such that

E\[Y \] = &beta;<sub>o</sub> +  &beta;<sub>X</sub> X<sub>true</sub> + &beta;<sub>G</sub> G

All of these values are inputted by the user (i.e. the intercept beta0, the effect of the exposure X<sub>true</sub> on the outcome as  &beta;<sub>X</sub>, and the effect of the SNP G directly on the outcome as  &beta;<sub>G</sub>).

After the SNP G, exposure X, and outcome Y are generated, then the reverseDirection function runs the MR Steiger approach to determine if the measured exposure X causes the outcome Y.

## Output
This function outputs the percent of simulations where the correct direction is detected between the exposure X and outcome Y using the MR Steiger approach. The R functions outputs the percent of simulations where the 3 cases detailed in (Hemani et al., 2017) are detected:

#### case 1: X->Y if the p-value from the Steiger correlation is less than alpha and p-value from the MR approach is less than alpha and the Steiger correlation Z>0
#### case 2: X<-Y if the p-value from the Steiger correlation is less than alpha and p-value from the MR approach is less than alpha and the Steiger correlation Z<0
#### case 3: inconclusive if the p-value from the Steiger correlation is greater than alpha or the p-value from the MR approach is greater than alpha 

The percent of simulations where the p-value from the Steiger correlation and MR are less than alpha are outputted (Steiger and MR, respectively). The correlation between the SNP G and the exposure X (corGX), correlation between the SNP G and the outcome Y (corGY), and the correlation between the exposure X and the outcome Y (corXY) are given.

## Example:
Consider an example with 100 subjects (input n=100) with a MAF of 50 (input MAF=0.5). Consider no pleiotropy or measurement error (input measurementError = F, pleiotropy = F). Then, let the exposure X be generated from a normal distribution with a variance of 1 (input varX = 1) and mean such that 
E\[X \] = 0 + 0.4 G
(input gamma0=0, gammaX=0.4). The outcome Y is generated from a normal distribution with a variance of 1 (input varY = 1) and mean such that 
E\[Y \] = 0 + &beta;<sub>X</sub> X 
(input beta0 = 0) and &beta;<sub>X</sub> varies from 0 to 1 by 0.25 (betaX = seq(from = 0, to = 1, length.out = 4)). The R code to run this example is given below.

```
library(reverseDirection)

results<-reverseDirection(nSim = 10000, n = 100, MAF = 0.5, gamma0 = 0, gammaG = 1, varX = 1, 
measurementError = F, delta0 = 0, deltaX = 1, varME = 0.2, 
beta0 = 0, betaX = c(seq(from = 0, to = 0.5, by=0.1),seq(from = 0.75, to = 2, by=0.25)), 
pleiotropy = F, betaG = 1, varY = 0.2, 
sig.level = 0.05, SEED = 1, plot.pdf = T, plot.name = "NoPleiotropyNoMeasurementError")


round(results$matrix,2)
```

The function outputs the following matrix and plot where each row corresponds to &beta;<sub>X</sub> (input betaX). As seen below, mostly case 3 is detected, which means that the MR Steiger method is inconclusive as to whether the exposure X causes the outcome Y.
```
      case1 case2 case3   Z+ Steiger   MR corGX corGY corXY
 [1,]  0.05     0  0.95 1.00    0.99 0.05  0.58  0.00  0.00
 [2,]  0.31     0  0.69 1.00    0.98 0.33  0.58  0.15  0.26
 [3,]  0.74     0  0.26 1.00    0.94 0.81  0.58  0.28  0.48
 [4,]  0.80     0  0.20 1.00    0.83 0.97  0.58  0.37  0.63
 [5,]  0.70     0  0.30 1.00    0.71 1.00  0.58  0.42  0.74
 [6,]  0.57     0  0.43 0.99    0.58 1.00  0.58  0.46  0.81
 [7,]  0.35     0  0.65 0.94    0.35 1.00  0.58  0.52  0.90
 [8,]  0.23     0  0.77 0.89    0.23 1.00  0.58  0.54  0.94
 [9,]  0.16     0  0.83 0.83    0.17 1.00  0.58  0.55  0.96
[10,]  0.12     0  0.87 0.79    0.13 1.00  0.58  0.56  0.97
[11,]  0.11     0  0.89 0.76    0.11 1.00  0.58  0.57  0.98
[12,]  0.09     0  0.90 0.74    0.10 1.00  0.58  0.57  0.98
```

<img src="reverseDirection.png" width="500">

## Reference:
Hemani G, Tilling K, Davey Smith G (2017) Orienting the causal relationship between imprecisely measured traits using GWAS summary data. PLOS Genetics 13(11): e1007081.
