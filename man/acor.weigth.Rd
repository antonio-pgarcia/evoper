% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abm-acor.R
\name{acor.weigth}
\alias{acor.weigth}
\title{Weight calculation for ant colony optimization}
\usage{
acor.weigth(q, k, l)
}
\arguments{
\item{q}{The Algorithm parameter. When small best-ranked solution is preferred}

\item{k}{The Archive size}

\item{l}{The lth element of algorithm solution archive T}
}
\value{
A scalar or a vector with calculated weigth.
}
\description{
Calculates the weight element of ACOr algorithm for the
solution archive.
}
\references{
[1] Socha, K., & Dorigo, M. (2008). Ant colony optimization for continuous domains.
European Journal of Operational Research, 185(3), 1155-1173.
http://doi.org/10.1016/j.ejor.2006.06.046
}
