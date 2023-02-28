#' STORM base columns
#'
#' Used internally by brainSTORM but may be useful
#'
#' @format character
#'
"storm_baseCols"

#' STORM base coordinate columns
#'
#' Used internally by brainSTORM but may be useful
#'
#' @format character
#'
"storm_baseCoorCols"

#' STORM base coordinate columns
#'
#' Used internally by brainSTORM but may be useful
#'
#' @format character
#'
"RNAmods_vec"

#' Pseudouridine scores
#'
#' @format character
#'
"Y_metrics"

#' 2-O-methylation scores
#'
#' @format character
#' @aliases Nm_scores
"Nm_metrics"

#' ac4C scores
#'
#' @format character
#' @aliases ac4C_scores
#'
"ac4C_metrics"

#' m1A scores
#' @aliases m1A_scores
#'
#' @format character
#'
"m1A_metrics"

#' m7G scores
#'
#' @format character
#' @aliases m7G_scores
#'
"m7G_metrics"

#' m5C scores
#'
#' @format character
#' @aliases m5C_scores
#'
"m5C_metrics"

#' m3U scores
#'
#' @format character
#' @aliases m3U_scores
#'
"m3U_metrics"

#' m1acp3Y metrics
#'
#' @format character
#'
"m1acp3Y_metrics"

#' m62A metrics
#'
#' @format character
#'
"m62A_metrics"

#' Yeast STORM object
#'
#' Example STORM object derived from STORMseq data of yeast.
#'
#' @format list
#'
"yeast_STORM"

#' Yeast rRNA modifications reference
#'
#' rRNA modifications in yeast according to Taoka et al.,
#' 2016 (Nucleic acids research)
#'
#' @format data.frame
#'
"rRNAmods_Sc_Taoka"

#' Variables list for experimental design table creation
#'
#' @format list
#'
"vList"

#' RNAmod metrics table
#'
#' Table with all RNAmod metrics used in brainSTORM with their respective RNAmod,
#' groups, RTase, and function used to generate it.
#'
#' @format data.frame
#'
"RNAmod_metrics"

#' All metrics vector
#'
#' A vector containing all brainSTORM metrics
#'
#' @format character
#'
"allMetrics"
