% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gLRTH_L.R
\name{gLRTH_L}
\alias{gLRTH_L}
\title{The function for the likelihood ratio test for genetic linkage under transmission heterogeneity}
\usage{
gLRTH_L(n0, n1, n2)
}
\arguments{
\item{n0}{Number of affected sibling pairs both of which inherited A from their heterozygous parent Aa}

\item{n1}{Number of affected sibling pairs which one inherited A and the other inherited a from their heterozygous parent Aa}

\item{n2}{Number of affected sibling pairs both of which inherited a from their heterozygous parent Aa}
}
\value{
The test statistic and asymptotic p-value for the likelihood ratio test for linkage analysis under genetic heterogeneity
}
\description{
We consider a binary trait and focus on detecting a transmission heterogeneity at a single locus with two
alleles \eqn{A} and \eqn{a}. We consider independent families each with one marker homozygous (\eqn{AA}) parent, one marker
heterozygous parent (\eqn{Aa}) and two diseased children. This likelihood ratio test is to test transmission
heterogeneity of preferential transmission of marker allele "a" to an affected child based on a binomial
mixture model with \eqn{J} components (\eqn{J \ge 2}),
\deqn{P_{\eta}(X_D=g)=\sum_{j=1}^J \alpha_j B_2(g, \theta_j), \; g=0, 1, 2, \; J \geq 2, \; \sum_{j=1}^J \alpha_j=1, \; \theta_j, \alpha_j \in (0, 1),}
where \eqn{\eta=(\eta_j)_{j \leq J}, \eta_j=(\theta_j, \alpha_j)^T, j=1, \ldots, J},
\eqn{B_2(g, \theta_j)} is the probability mass function for a binomial distribution \eqn{X \sim Bin(2, \theta_j)},
and \eqn{\theta_i=\theta_j} if and only if \eqn{i=j}. \eqn{\theta_j} is the probability
of transmission of the allele of interest in a subgroup of families \eqn{j}.
In particular, \eqn{J} is likely to be quite
large for many of the complex disease under transmission heterogeneity. Note that this LRT can
be applied to genome-wide linkage analysis without the need to know the exact value of \eqn{J} while allowing \eqn{J \ge 2}.
}
\examples{
gLRTH_L(n0=100, n1=70, n2=30)
}
\references{
Shao Y. (2014) Linkage analysis, originally published
in Encyclopedia of Quantitative Risk Analysis and Assessment, John Wiley & Sons, Ltd, USA, 2008, and
republished in Wiley StatsRef: Statistics Reference Online 2014.
}
\author{
Xiaoxia Han and Yongzhao Shao
}
