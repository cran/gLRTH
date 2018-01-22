#' The function for the likelihood ratio test for genome-wide association under genetic
#' heterogeneity with genotype frequencies as input values
#'
#' @description  We consider a binary trait and focus on detecting association with disease at a single locus with two
#' alleles \eqn{A} and \eqn{a}. The likelihood ratio test is based
#' on a binomial mixture model of \eqn{J} components (\eqn{J \ge 2}) for diseased cases:
#' \deqn{P_{\eta}(X_D=g)=\sum_{j=1}^J \alpha_j B_2(g, \theta_j), \; g=0, 1, 2, \; J \geq 2, \; \sum_{j=1}^J \alpha_j=1, \; \theta_j, \alpha_j \in (0, 1),}
#' where \eqn{\eta=(\eta_j)_{j \leq J}, \eta_j=(\theta_j, \alpha_j)^T, j=1, \ldots, J},
#' \eqn{B_2(g, \theta_j)} is the probability mass function for a binomial distribution \eqn{X \sim Bin(2, \theta_j)},
#' and \eqn{\theta_i=\theta_j} if and only if \eqn{i=j}. \eqn{\theta_j} is the probability
#' of having the allele of interest on one chromosome for a subgroup of case \eqn{j}.
#' In particular, \eqn{J} is likely to be quite
#' large for many of the complex disease with genetic heterogeneity. Note that the LRT-H can
#' be applied to association studies without the need to know the exact value of \eqn{J} while allowing \eqn{J \ge 2}.
#' @return The test statistic and asymptotic p-value for the likelihood ratio test for GWAS under genetic heterogeneity
#' @param n0 AA genotype frequency in case
#' @param n1 Aa genotype frequency in case
#' @param n2 aa genotype frequency in case
#' @param m0 AA genotype frequency in control
#' @param m1 Aa genotype frequency in control
#' @param m2 aa genotype frequency in control
#' @author Xiaoxia Han and Yongzhao Shao
#' @references
#' Qian M., Shao Y. (2013) A Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity.
#' Annals of Human Genetics, 77(2): 174-182.
#' @examples
#' gLRTH_A(n0=2940, n1=738, n2=53, m0=3601, m1=1173, m2=117)
#' @export
#' @importFrom stats pchisq

gLRTH_A<-function(n0, n1, n2, m0, m1, m2){
  n<-n0+n1+n2
  m<-m0+m1+m2
  p0<-(n2+m2+(n1+m1)/2)/(n+m)
  pH<-(m2+m1/2)/m
  pD<-(n2+n1/2)/n

  LRT.stat<-ifelse(4*n0*n2>n1^2,
                   2*log(((1-pH)/(1-p0))^(2*m0)*(pH*(1-pH)/p0/(1-p0))^m1*(pH/p0)^(2*m2)*
                           (n0/(n*(1-p0)^2))^n0*(n1/(n*2*p0*(1-p0)))^n1*(n2/(n*p0^2))^n2),
                   2*log(((1-pH)/(1-p0))^(2*m0)*(pH*(1-pH)/p0/(1-p0))^m1*(pH/p0)^(2*m2)*
                           ((1-pD)/(1-p0))^(2*n0)*(pD*(1-pD)/p0/(1-p0))^n1*(pD/p0)^(2*n2)))

  pval<-(pchisq(LRT.stat, 1, lower.tail = F) + pchisq(LRT.stat, 2, lower.tail = F))/2

  return(list(chisq.stat=LRT.stat, pval=pval))
}
