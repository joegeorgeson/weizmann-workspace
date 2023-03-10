% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{add_ER}
\alias{add_ER}
\title{Add read-end to coverage ratio}
\usage{
add_ER(
  STORM,
  group_A,
  newColName = "auto",
  onNucs = c("A", "C", "G", "T"),
  minCov = 50,
  nCores = 1
)
}
\arguments{
\item{STORM}{list. STORM object as output by \code{\link{storm_STORM}}}

\item{group_A}{character. Name of group to be compared as found in STORM$META$group}

\item{newColName}{character. Name of calculated metric to be stored in STORM$RES,
assigned by default based on group_A and group_B}

\item{onNucs}{character. Nucleotide(s) in which the metric will be calculated}

\item{minCov}{integer. Minimal coverage in position for metric to be calculated}
}
\description{
Add read-end to coverage ratio
}
