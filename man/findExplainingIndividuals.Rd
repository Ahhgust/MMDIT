% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{findExplainingIndividuals}
\alias{findExplainingIndividuals}
\title{Find sets of individuals that can explain the graph.}
\usage{
findExplainingIndividuals(sgrap, toExS, numInds, maxMissNode)
}
\arguments{
\item{sgrap}{The graph.}

\item{toExS}{The graph traversals.}

\item{numInds}{The number of individuals to consider.}

\item{maxMissNode}{The maximum number of missed nodes.}
}
\value{
The sets of individuals that can explain the graph.
}
\description{
Find sets of individuals that can explain the graph.
}
