% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{add_CytPer}
\alias{add_CytPer}
\title{Add cytidine persistance index to STORM object}
\usage{
add_CytPer(
  STORM,
  group_A,
  group_B,
  newColName = "auto",
  onNucs = c("C"),
  minReads = 30,
  nCores = 1
)
}
\arguments{
\item{STORM}{list. STORM object as output by \code{\link{storm_STORM}}}

\item{group_A}{character. Name of group to be compared as found in STORM$META$group}

\item{group_B}{character. Name of group to be compared as found in STORM$META$group}

\item{newColName}{character. Name of calculated metric to be stored in STORM$RES,
assigned by default based on group_A and group_B}

\item{onNucs}{character. Nucleotide(s) in which the metric will be calculated}

\item{minReads}{numeric. Minimum number of nucleotides read in position to
compute the output value.}
}
\description{
Cytidine score modified from (Squires et al. 2012). It accounts for positions
with Thymine reads in Cytidine positions in the Control/Mock sample.
CytPer = (C_bisulphite / (C_bisulphite + T_bisulphite)) - (1 - (C_control / C_control + T_control))
}
