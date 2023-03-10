% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{add_DR}
\alias{add_DR}
\title{Add deletion rate}
\usage{
add_DR(
  STORM,
  group_A,
  newColName = "auto",
  onNucs = c("A", "C", "G", "T"),
  minReads = 30,
  nCores = 1
)
}
\arguments{
\item{STORM}{list. STORM object as output by \code{\link{storm_STORM}}}

\item{group_A}{character. Name of group to be compared as found in STORM$META$group}

\item{newColName}{character. Name of calculated metric to be stored in STORM$RES,
assigned by default based on group_A and group_B}

\item{onNucs}{character. Nucleotide(s) in which the metric will be calculated}

\item{minReads}{numeric. Minimum number of nucleotides read in position to
compute the output value.}
}
\description{
Adds the deletion to nucleotide reads ratio
}
