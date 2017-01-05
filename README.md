An R package for the Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity
================

Example
-------

``` r
library(gLRTH)
# input is 3 by 2 genotype summary table
# n0, n1 and n2 are AA, Aa and aa frequences in case
# m0, m1, and m2 are Aa Aa and aa frequences in control
gLRTH(n0=2940, n1=738, n2=53, m0=3601, m1=1173, m2=117)

# input is genotype and disease status of individual level 
# simulate genotype for 200 case and 200 control
disease<-c(rep(1, 200), rep(0, 200))
geno1<-c(rbinom(n=50, size=2, prob=0.5), rbinom(n=150, size=2, prob=0.23))
geno2<-rbinom(n=200, size=2, prob=0.5)
geno<-c(geno1, geno2)
gLRTH2(g=geno, d=disease)
```
