% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abm-ees1.R
\name{ees1.selection}
\alias{ees1.selection}
\title{ees.selection}
\usage{
ees1.selection(s0, s1, kkappa)
}
\arguments{
\item{s0}{The current best solution set}

\item{s1}{The new solution}

\item{kkappa}{The selection pressure}
}
\description{
Select the elements with best fitness but accept uphill
moves with probability 'kkappa'.
}
