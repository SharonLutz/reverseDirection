# reverseDirection
Examines the MR Steiger approach to detect the directionality between the mediator and outcome through simulation studies

## Installation
```
# you need R v3.4 or higher
install.packages("devtools")  # devtools must be installed first
install.packages("MendelianRandomization") 
install.packages("psych")
devtools::install_github("MRCIEU/TwoSampleMR") 

devtools::install_github("SharonLutz/reverseDirection")
```

## Input
nSNP is the number of SNPs generated from a binomial distribution for n subjects (input n) for a given minor allele frequency (input vector MAF).

For the SNPs Xi, the mediator/ exposure M is generated from a normal distribution with the variance (input varM) and the mean as follows:

E\[M\] = &gamma;<sub>o</sub> + &sum; &gamma;<sub>X</sub>  X<sub>i</sub> 

All of these values are inputted by the user (i.e. the intercept gamma0, and the genetic effect size as a vector gammaX).

The outcome Y is generated from a normal distribution with the variance (input varY) and the mean as follows:

E\[Y\] = &beta;<sub>o</sub> +  &beta;<sub>M</sub> M 

All of these values are inputted by the user (i.e. the intercept beta0 and the effect of the mediator directly on the outcome as betaM).

After the SNPs X, mediator M, and outcome Y are generated, then the reverseDirection function runs the MR Steiger approach to determine if the mediator M causes the outcome Y.


## Example:

```
library(reverseDirection)
?reverseDirection # For details on this function

reverseDirection(n = 1000)

```

## Output:
<img src="reverseDirectionplot.png" width="500">


