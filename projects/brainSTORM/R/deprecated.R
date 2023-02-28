#' @export
storm_add_scoreA5p <- function(STORM, group_A, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    lifecycle::deprecate_warn("0.0.6", "storm_add_scoreA5p()", "add_scoreA5p()")
    add_scoreA5p(STORM, group_A, newColName, onNucs, minMedCov)
}

#' @export
storm_add_scoreA3p <- function(STORM, group_A, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    lifecycle::deprecate_warn("0.0.6", "storm_add_scoreA3p()", "add_scoreA3p()")
    add_scoreA3p(STORM, group_A, newColName, onNucs, minMedCov)
}

#' @export
storm_makeCalls <- function(STORM){
    lifecycle::deprecate_warn("0.0.6", "storm_makeCalls()", "storm_assignScores()")
    storm_assignScores(STORM)
}


#' Add all default metrics for STORMseq v1
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated as the STORMseq library treatments were expanded.
#' This functions applied only for the library treatments of lib 499.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#'
#' @return
#'
#' @examples
add_default_metrics_v1 <- function(STORM){
    lifecycle::deprecate_warn("0.0.6", "add_default_metrics_v1()", "add_default_metrics()")
        STORM %>%
        add_SRD1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRD1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_SRlog2FCh1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRlog2FCh1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_scoreA3p("Mock.TGIRT") %>%
        add_scoreA3p("Mock.SSIII") %>%
        add_2OmeScore("MocklowdNTPs.SSIII", "Mock.SSIII", perNuc_Zscaling = T) %>%
        add_CtoT_MRD("Ac4C.TGIRT", "Mock.TGIRT") %>%
        add_CtoT_MRD("Ac4C.TGIRT", "DeacetylatedAc4C.TGIRT") %>%
        add_CtoT_MRD("Ac4C.SSIII", "Mock.SSIII") %>%
        add_CtoT_MRD("Ac4C.SSIII", "DeacetylatedAc4C.SSIII") %>%
        add_SRD1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_SRlog2FCh1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_MRD("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRD1bpDS("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRlog2FCh1bpDS("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRD1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRD1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_MRD("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_CytPer("m5C.TGIRT", "Mock.TGIRT") %>%
        add_MRD("Mock.TGIRT", "AlkBmix.TGIRT")
}
