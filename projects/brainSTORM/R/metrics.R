# PREDEFINED METRICS ###########################################################

#' Add all default metrics for STORMseq
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param nCores
#'
#' @return
#' @export
#'
add_default_metrics <- function(STORM, nCores = 1){
    STORM %>%
        add_2OmeScore("MocklowdNTPs.SSIII", "Mock.SSIII", nCores = nCores) %>%
        add_CtoT_MRD("Ac4C.SSIII", "DeacetylatedAc4C.SSIII", nCores = nCores) %>%
        add_CtoT_MRD("Ac4C.SSIII", "Mock.SSIII", nCores = nCores) %>%
        add_CtoT_MRD("Ac4C.TGIRT", "DeacetylatedAc4C.TGIRT", nCores = nCores) %>%
        add_CtoT_MRD("Ac4C.TGIRT", "Mock.TGIRT", nCores = nCores) %>%
        add_CytPer("m5C.TGIRT", "Mock.TGIRT", nCores = nCores) %>%
        add_CytPer("RBSseqHeatMg.SSIII", "Mock.SSIII", nCores = nCores) %>%
        add_CytPer("RBSseqHeatMg.TGIRT", "Mock.TGIRT", nCores = nCores) %>%
        add_DRD("m5C.TGIRT", "Mock.TGIRT", nCores = nCores) %>%
        add_DRD("RBSseqHeatMg.SSIII", "Mock.SSIII", nCores = nCores) %>%
        add_DRD("RBSseqHeatMg.TGIRT", "Mock.TGIRT", nCores = nCores) %>%
        add_MRD("DeacetylatedAc4C.SSIII", "Ac4C.SSIII", nCores = nCores) %>%
        add_MRD("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT", nCores = nCores) %>%
        add_MRD("Mock.SSIII", "RBSseqHeatMg.SSIII", nCores = nCores) %>%
        add_MRD("Mock.TGIRT", "AlkBmix.TGIRT", nCores = nCores) %>%
        add_MRD("Mock.TGIRT", "Dimroth.TGIRT", nCores = nCores) %>%
        add_MRD("Mock.TGIRT", "RBSseqHeatMg.TGIRT", nCores = nCores) %>%
        add_MRD("NaBH4HydBiotin.RTHIV", "Mock.RTHIV", nCores = nCores) %>%
        add_MRD("NaBH4HydBiotin.TGIRT", "Mock.TGIRT", nCores = nCores) %>%
        add_MRDns(group_A = "Mock.TGIRT", group_B = "m5C.TGIRT", refNuc = "A", misNuc = "G", nCores = nCores) %>%
        add_MRDns(group_A = "Mock.SSIII", group_B = "m5C.SSIII", refNuc = "A", misNuc = "G", nCores = nCores) %>%
        add_scoreA3p("Mock.RTHIV", nCores = nCores) %>%
        add_scoreA3p("Mock.SSIII", nCores = nCores) %>%
        add_scoreA3p("Mock.TGIRT", nCores = nCores) %>%
        add_scoreC3p("Mock.RTHIV", nCores = nCores) %>%
        add_scoreC3p("Mock.SSIII", nCores = nCores) %>%
        add_scoreC3p("Mock.TGIRT", nCores = nCores) %>%
        add_SRD1bpDS("CMC.SSIII", "Mock.SSIII", nCores = nCores) %>%
        add_SRD1bpDS("CMC.TGIRT", "Mock.TGIRT", nCores = nCores) %>%
        add_SRD1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII", nCores = nCores) %>%
        add_SRD1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT", nCores = nCores) %>%
        add_SRD1bpDS("Mock.SSIII", "Dimroth.SSIII", nCores = nCores) %>%
        add_SRlog2FCh1bpDS("CMC.SSIII", "Mock.SSIII", nCores = nCores) %>%
        add_SRlog2FCh1bpDS("CMC.TGIRT", "Mock.TGIRT", nCores = nCores) %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII", nCores = nCores) %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT", nCores = nCores) %>%
        add_SRlog2FCh1bpDS("Mock.SSIII", "Dimroth.SSIII")
}

add_Y_metrics <- function(STORM){
    STORM %>%
        add_SRD1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRD1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_SRlog2FCh1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRlog2FCh1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_DRD("RBSseqHeatMg.TGIRT", "Mock.TGIRT") %>%
        add_DRD("RBSseqHeatMg.SSIII", "Mock.SSIII") %>%
        add_DRD("m5C.TGIRT", "Mock.TGIRT")
}

add_Nm_metrics <- function(STORM){
    STORM %>%
        add_scoreA3p("Mock.TGIRT") %>%
        add_scoreA3p("Mock.SSIII") %>%
        add_scoreA3p("Mock.RTHIV") %>%
        add_2OmeScore("MocklowdNTPs.SSIII", "Mock.SSIII", perNuc_Zscaling = T)
}

add_ac4C_metrics <- function(STORM){
    STORM %>%
        add_CtoT_MRD("Ac4C.TGIRT", "Mock.TGIRT") %>%
        add_CtoT_MRD("Ac4C.SSIII", "Mock.SSIII") %>%
        add_CtoT_MRD("Ac4C.TGIRT", "DeacetylatedAc4C.TGIRT") %>%
        add_CtoT_MRD("Ac4C.SSIII", "DeacetylatedAc4C.SSIII")
}

add_m1A_metrics <- function(STORM){
    STORM %>%
        add_SRD1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_SRlog2FCh1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_MRD("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_MRD("Mock.TGIRT", "RBSseqHeatMg.TGIRT") %>%
        add_MRD("Mock.SSIII", "RBSseqHeatMg.SSIII")
}

add_m7G_metrics <- function(STORM){
    STORM %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRD1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRD1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_MRD("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("NaBH4HydBiotin.TGIRT", "Mock.TGIRT") %>%
        add_MRD("NaBH4HydBiotin.RTHIV", "Mock.RTHIV")
}

add_m5C_metrics <- function(STORM){
    STORM %>%
        add_CytPer("m5C.TGIRT", "Mock.TGIRT") %>%
        add_CytPer("RBSseqHeatMg.TGIRT", "Mock.TGIRT") %>%
        add_CytPer("RBSseqHeatMg.SSIII", "Mock.SSIII")
}

add_m3U_metrics <- function(STORM){
    STORM %>%
        add_MRD("Mock.TGIRT", "AlkBmix.TGIRT")
}

add_Nm_MgR_metrics <- function(STORM){
    STORM %>%
        add_ER1bpDS("MgR") %>%
        add_ER1bpDS("MgR-OED") %>%
        add_ER1bpDS("Mock") %>%
        add_scoreZ_3p("MgR") %>%
        add_scoreZ_3p("MgR-OED") %>%
        add_scoreZ_3p("Mock") %>%
        add_ERlog2FCh("MgR", "Mock") %>%
        add_ERlog2FCh1bpDS("MgR", "Mock") %>%
        add_ERlog2FCh("MgR-OED", "Mock") %>%
        add_ERlog2FCh1bpDS("MgR-OED", "Mock")
}

add_m62A_metrics <- function(STORM){
    STORM %>%
        add_MRDns(group_A = "Mock.TGIRT", group_B = "m5C.TGIRT", refNuc = "A", misNuc = "G") %>%
        add_MRDns(group_A = "Mock.SSIII", group_B = "m5C.SSIII", refNuc = "A", misNuc = "G")
}


# METRIC ASSIGNMENT FUNCTIONS ##################################################

#' Initialize STORM$CALLS
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#'
#' @return
#' @export
#'
#' @examples
hlp_start_CALLS <- function(STORM){
    RNAmods <- check_whichRNAmods(STORM)
    STORM$CALLS <- lapply(seq_along(unique(STORM$META$set)), function(x) {
        lapply(RNAmods, function(x) NULL) %>% magrittr::set_names(RNAmods)
    })
    names(STORM$CALLS) <- unique(STORM$META$set)
    STORM
}

#' Assign scores to respective RNAmods in STORM$CALLS
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}.
#' @param RNAmod character. Name of RNAmod to assign scores to.
#' @param scores character. Names of scores to be assigned to RNAmod.
#' @param onNucs character. Vector of base nucleotides to include in scores, all
#' nucleotides by default, i.e. c("A", "C", "G", "T").
#'
#' @return
#' @export
#'
#' @examples
hlp_assign_scores <- function(STORM, RNAmod, scores, onNucs = c("A", "C", "G", "T")){
    if(any(scores %in% unique(STORM$RES$metric))){
        tmp <- STORM$RES %>%
            tidyr::pivot_wider(names_from = metric, values_from = score) %>%
            data.table::data.table()
        scoresInTable <- intersect(scores, colnames(tmp))
        selVars <- c(storm_baseCols, scoresInTable)
        tmp <- tmp[,names(tmp) %in% selVars, with = FALSE]
        tmp[!(tmp$refSeq %in% onNucs), scoresInTable] <- NA
        scoresNotInTable <- setdiff(scores, colnames(tmp))
        if(length(scoresNotInTable) > 0){
            warning(paste0(scoresNotInTable, collapse = " "), " metric(s) are not ",
                    "present in STORM$RES$metric, and were not added to CALLS. ",
                    "This is just a warning.")
        }
        tmp <- split(tmp, tmp$set)
        for(iN in names(tmp)){
            STORM$CALLS[[iN]][[RNAmod]] <- tmp[[iN]]
        }
        return(STORM)
    }else{
        warning(RNAmod, " was not added in CALLS as no relevant metric is present.")
        STORM
    }
}


#' Assigning default Scores to tables by sets of samples
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#'
#' @return
#' @export
#'
#' @examples
storm_assignScores <- function(STORM){
    hlp_start_CALLS(STORM) %>%
        hlp_assign_scores("Y", Y_scores) %>%
        hlp_assign_scores("Am", Nm_scores, "A") %>%
        hlp_assign_scores("Tm", Nm_scores, "T") %>%
        hlp_assign_scores("Gm", Nm_scores, "G") %>%
        hlp_assign_scores("Cm", Nm_scores, "C") %>%
        hlp_assign_scores("m5C", m5C_scores) %>%
        hlp_assign_scores("ac4C", ac4C_scores) %>%
        hlp_assign_scores("m1A", m1A_scores) %>%
        hlp_assign_scores("m7G", m7G_scores) %>%
        hlp_assign_scores("m3U", m3U_scores) %>%
        hlp_assign_scores("m1acp3Y", m1acp3Y_metrics) %>%
        hlp_assign_scores("m62A", m62A_metrics)
}
# BASE FUNCTIONS ###############################################################

# Between two groups ###########################################################

#' Add End rate log2-Fold-Change 1bp downstream
#'
#' Calculates the log2-Fold-Change of the end rate ratio 1bp downstream
#' of group_A over that of group_B per set of samples in a STORM object.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
add_ERlog2FCh1bpDS=function(STORM, group_A, group_B, newColName = "auto",
                            onNucs = c("A", "C", "G", "T"), minCov = 50, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("ERlog2FCh1bpDS", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- add_EndRate1bpDS(STORM$DATA[[id_A]], minCov = minCov)
        DT_B <- add_EndRate1bpDS(STORM$DATA[[id_B]], minCov = minCov)
        tmpRES <- data.table::data.table(log2(DT_A$endRate_1bpDS / DT_B$endRate_1bpDS))
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add End rate log2-Fold-Change
#'
#' Calculates the log2-Fold-Change of the end rate ratio
#' of group_A over that of group_B per set of samples in a STORM object.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
add_ERlog2FCh <- function(STORM, group_A, group_B, newColName = "auto",
                          onNucs = c("A", "C", "G", "T"), minCov = 50, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("ERlog2FCh", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- add_EndRate(STORM$DATA[[id_A]], minCov = minCov)
        DT_B <- add_EndRate(STORM$DATA[[id_B]], minCov = minCov)
        tmpRES <- data.table::data.table(log2(DT_A$end_Ratio / DT_B$end_Ratio))
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}


#' Add Start Rate Difference
#'
#' Calculates the start ratio difference of group_A minus that
#' of group_B per set of samples in a STORM object.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
add_SRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T"), minCov = 50, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("SRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- add_StartRate(STORM$DATA[[id_A]], minCov = minCov)
        DT_B <- add_StartRate(STORM$DATA[[id_B]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$startRatio - DT_B$startRatio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add Start Rate Difference 1bp down-stream
#'
#' Calculates the start ratio difference 1bp down-stream of group_A minus that
#' of group_B per set of samples in a STORM object.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
add_SRD1bpDS <- function(STORM, group_A, group_B, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minCov = 50, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("SRD1bpDS", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]], minCov = minCov)
        DT_B <- add_StartRate_1bpDS(STORM$DATA[[id_B]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$startRate_1bpDS - DT_B$startRate_1bpDS)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add Start ratio log2-Fold-Change 1bp down-stream
#'
#' Calculates the log2-Fold-Change of the start ratio difference 1bp down-stream
#' of group_A over that of group_B per set of samples in a STORM object.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_SRlog2FCh1bpDS <- function(STORM, group_A, group_B, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minCov = 50,
                               nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("SRlog2FCh1bpDS", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]], minCov = minCov)
        DT_B <- add_StartRate_1bpDS(STORM$DATA[[id_B]], minCov = minCov)
        tmpRES <- data.table::data.table(log2(DT_A$startRate_1bpDS / DT_B$startRate_1bpDS))
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add 2-O-methylation score to STORM-seq object
#'
#' Methylation score that considers the log2-FCH of start ratio between
#' two RNA-seq libraries (low-dNTPs and high-dNTPs), subtracting the "noise"
#' from a surrounding window.
#'
#' Incarnato et al. 2017: For Nm Incarnato and collabs. used the fold change
#' as the log2 of the ratio between the stoppage ratio at the low dNTPs
#' concentration and the stoppage ratio in the high dNTPs concentration sample.
#'  Because the stoppage ratio seems to be strongly region specific, they
#'  subtracted the local background defined as the mean ratio of the stoppage
#'  in the neighborhood of a given position (+/- 5 nucleotides, excluding
#'  given position).
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param lib_LowdNTPs character. Name of group with low-dNTPs library preparation
#' found in STORM$META$group
#' @param lib_HighdNTPs character. Name of group with high-dNTPs library preparation
#' found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param flankSize integer. Length of window to each side of position to
#' consider for metric calculation
#' @param perNuc_Zscaling
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_2OmeScore <- function(STORM,
                          lib_LowdNTPs,
                          lib_HighdNTPs,
                          newColName = "auto",
                          flankSize = 5,
                          perNuc_Zscaling = F,
                          onNucs = c("A", "C", "G", "T"),
                          minCov = 50,
                          nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("NmStopScore", lib_LowdNTPs, lib_HighdNTPs, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == lib_LowdNTPs,]$set,
                      STORM$META[STORM$META$group == lib_HighdNTPs,]$set)
    if(length(sets) == 0){warning("No sets found with both ", lib_LowdNTPs,
                                  " and ", lib_HighdNTPs, " group labels.\n",
                                  newColName, " was not calculated.")
        return(STORM)}
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == lib_LowdNTPs,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == lib_HighdNTPs,]$id
        DT <- detect_2OmeScore(DT_low = STORM$DATA[[id_A]],
                               DT_high = STORM$DATA[[id_B]],
                               neighFlankSize = flankSize,
                               minCov = minCov)
        tmpRES <- data.table::data.table(twoOme_score = DT$twoOme_score)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT[,c("chr", "gencoor", "strand", "gene",
                                               "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    })
    if(perNuc_Zscaling == T){
        OUT <- lapply(OUT, function(DT){
            iNucs <- DT$refSeq %>% unique()
            for(iN in iNucs){
                selPos <- which(DT$refSeq == iN)
                DT$score[selPos] <- scale(DT$score[selPos])
            }
            DT
        })
    }
    OUT <- do.call(OUT, what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add Misincorporation Rate Difference
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#'
#' @return
#' @export
#'
#' @examples
add_MRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T"), minReads = 30, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("MRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- STORM$DATA[[id_A]] %>% add_diffNucToRefRatio(minReads = minReads)
        DT_B <- STORM$DATA[[id_B]] %>% add_diffNucToRefRatio(minReads = minReads)
        tmpRES <- data.table::data.table(DT_A$diffToRef_Ratio - DT_B$diffToRef_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add Nucleotide Specific Misincorporation Rate Difference
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param refNuc character. Reference nucleotide
#' @param misNuc character. Misincorporated nucleotide
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#'
#' @return
#' @export
#'
#' @examples
add_MRDns <- function(STORM, group_A, group_B, refNuc, misNuc, minReads = 30, nCores = 1){
    newColName <- paste0("MRD", refNuc, "to", misNuc, "_", group_A, "_", group_B)
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        MR_colName <- paste0("MR_", refNuc, "to", misNuc)
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- STORM$DATA[[id_A]] %>%
            txtools::tx_add_misincRateNucSpec(
                refNuc = refNuc, misNuc = misNuc, minNucReads = minReads)
        DT_B <- STORM$DATA[[id_B]] %>%
            txtools::tx_add_misincRateNucSpec(
                refNuc = refNuc, misNuc = misNuc, minNucReads = minReads)
        tmpRES <- data.table::data.table(DT_A[[MR_colName]] - DT_B[[MR_colName]])
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% refNuc), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add cytidine persistance index to STORM object
#'
#' Cytidine score modified from (Squires et al. 2012). It accounts for positions
#' with Thymine reads in Cytidine positions in the Control/Mock sample.
#' CytPer = (C_bisulphite / (C_bisulphite + T_bisulphite)) - (1 - (C_control / C_control + T_control))
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#'
#' @return
#' @export
#'
#' @examples
add_CytPer <- function(STORM, group_A, group_B, newColName = "auto",
                       onNucs = c("C"), minReads = 30, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("CytPer", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- STORM$DATA[[id_A]] %>% detect_m5C(minReads = minReads)
        DT_B <- STORM$DATA[[id_B]] %>% detect_m5C(minReads = minReads)
        tmpRES <- data.table::data.table(DT_A$det_m5C - (1 - DT_B$det_m5C))
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add C-to-T misincorporation difference
#'
#' Calculates the difference of T conversion at C positions by subtracting the
#' result of "T"/("T"+"C") between groups. Commonly NaCNBH3-Treated vs
#' Control or vs Deacetylated negative control.
#'
#' The reaction of ac4C with sodium cyanoborohydride (NaCNBH3) under acidic
#' conditions forms the reduced nucleobase N4-acetyltetrahydrocytidine.
#' The altered structure of this reduced nucleobase compared with ac4C causes
#' the incorporation of non-cognate deoxynucleotide triphosphates (dNTPs)
#' upon reverse transcription (Sas-Chen et al. 2020).
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#'
#' @return
#' @export
#'
#' @examples
add_CtoT_MRD <- function(STORM, group_A, group_B, newColName = "auto",
                         onNucs = c("C"), minReads = 30, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("CtoT.MRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- STORM$DATA[[id_A]] %>% detect_ac4C(minReads = minReads)
        DT_B <- STORM$DATA[[id_B]] %>% detect_ac4C(minReads = minReads)
        tmpRES <- data.table::data.table(DT_A$det_ac4C - DT_B$det_ac4C)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add deletion rate difference
#'
#' Calculated the deletion rate to nucleotide reads for each position and
#' subtracts it for the rate in the control group.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#'
#' @return
#' @export
#'
#' @examples
add_DRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T"), minReads = 30, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("DRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set,
                      STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- STORM$DATA[[id_A]] %>% deletion_rate(minReads = minReads)
        DT_B <- STORM$DATA[[id_B]] %>% deletion_rate(minReads = minReads)
        tmpRES <- data.table::data.table(DT_A$deletion_Ratio - DT_B$deletion_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}



# One group metrics ############################################################

#' Add read-end rate 1bp downstream
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#'
#' @return
#' @export
#'
#' @examples
add_ER1bpDS=function(STORM, group_A, newColName = "auto",
                     onNucs = c("A", "C", "G", "T"), minCov = 50, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("ER1bpDS", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_EndRate1bpDS(STORM$DATA[[id_A]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$endRate_1bpDS)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add read-start to coverage rate
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_SR <- function(STORM, group_A, newColName = "auto",
                   onNucs = c("A", "C", "G", "T"), minCov = 50, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("SR", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_StartRate(STORM$DATA[[id_A]])
        tmpRES <- data.table::data.table(DT_A$start_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add read-start to coverage rate 1 bp down-stream
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_SR1bpDS <- function(STORM, group_A, newColName = "auto",
                        onNucs = c("A", "C", "G", "T"), minCov = 50, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("SR1bpDS", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$startRate_1bpDS)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add read-end to coverage ratio
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_ER <- function(STORM, group_A, newColName = "auto",
                   onNucs = c("A", "C", "G", "T"), minCov = 50, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("ER", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_EndRate(STORM$DATA[[id_A]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$end_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add ScoreA accounting for read-ends (3prime)
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param flankSize integer. Length of window to each side of position to
#' consider for metric calculation
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#'
#' @return
#' @export
#' @aliases storm_add_scoreA3p
#'
#' @examples
add_scoreA3p <- function(STORM, group_A, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minReads = 30,
                         flankSize = 6, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("ScoreA3p", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_scoreA_3p(STORM$DATA[[id_A]], flankSize = flankSize, minReads = minReads)
        tmpRES <- data.table::data.table(DT_A$scoreA_3p)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add ScoreA accounting for read-starts (5prime)
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param flankSize integer. Length of window to each side of position to
#' consider for metric calculation
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#'
#' @return
#' @export
#'
#' @examples
#' @aliases storm_add_scoreA5p
add_scoreA5p <- function(STORM, group_A, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minReads = 30,
                         flankSize = 6, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("ScoreA5p", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_scoreA_5p(STORM$DATA[[id_A]], minReads = minReads,
                              flankSize = flankSize)
        tmpRES <- data.table::data.table(DT_A$scoreA_5p)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[, c("chr", "gencoor", "strand", "gene",
                                                  "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add ScoreC accounting for read-stops (3prime)
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minMedCov integer. Minimal median coverage of window around position,
#' delimited by flankSize
#' @param flankSize integer. Length of window to each side of position to
#' consider for metric calculation
#' @param minReads integer. Minimum number of reads needed at site neighborhood
#' to calculate the score.
#'
#' @return
#' @export
#'
#' @examples
#' @aliases storm_add_scoreC3p
add_scoreC3p <- function(STORM, group_A, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minReads = 30,
                         flankSize = 6, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("ScoreC3p", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_scoreC_3p(STORM$DATA[[id_A]], minReads = minReads, flankSize = flankSize)
        tmpRES <- data.table::data.table(DT_A$scoreC)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[, c("chr", "gencoor", "strand", "gene",
                                                  "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add scoreZ accounting for read-ends (3prime)
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#' @param flankSize integer. Length of window to each side of position to
#' consider for metric calculation
#'
#' @return
#' @export
#'
#' @examples
add_scoreZ_3p = function(STORM, group_A, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minReads = 30,
                         flankSize = 6, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("scoreZ_3p", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_scoreZ_3p_util(STORM$DATA[[id_A]], minReads = minReads,
                                   flankSize = flankSize)
        tmpRES <- data.table::data.table(DT_A$scoreZ_3p)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add deletion rate
#'
#' Adds the deletion to nucleotide reads ratio
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#'
#' @return
#' @export
#'
#' @examples
add_DR <- function(STORM, group_A, newColName = "auto",
                   onNucs = c("A", "C", "G", "T"), minReads = 30, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("DR", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- deletion_rate(STORM$DATA[[id_A]], minReads = minReads)
        tmpRES <- data.table::data.table(DT_A$deletion_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add Misincorporation Rate
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minReads numeric. Minimum number of nucleotides read in position to
#' compute the output value.
#'
#' @return
#' @export
#'
#' @examples
add_MR <- function(STORM, group_A, newColName = "auto",
                   onNucs = c("A", "C", "G", "T"), minReads = 30, nCores = 1){
    if(newColName == "auto"){
        newColName <- paste("MR", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- parallel::mclapply(mc.cores = nCores, sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_diffNucToRefRatio(STORM$DATA[[id_A]], minReads = minReads)
        tmpRES <- data.table::data.table(DT_A$diffToRef_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Helper functions #############################################################
# Function for cytidine persistence to Bisulphite treatment
detect_m5C <- function(DT, minReads = 30){
    DT <- data.table::data.table(DT)
    tmpR <- (DT$`T` + DT$C)
    tmp <- round(DT$C / tmpR, 6)
    tmp[tmpR < minReads] <- NA
    tibble::add_column(DT, det_m5C = tmp)
}

# Function for cytidine persistence to Bisulphite treatment
detect_ac4C <- function(DT, minReads = 30){
    DT <- data.table::data.table(DT)
    tmpR <- (DT$`T` + DT$C)
    tmp <- round(DT$`T` / tmpR, 6)
    tmp[tmpR < minReads] <- NA
    tibble::add_column(DT, det_ac4C = tmp)
}

# Add column of Different Nucleotides to reference ratio, diffToRef and
# nucTotal columns are required
add_diffNucToRefRatio <- function(DT, minReads = 30){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("diffToRef_Ratio")
    DT <- txtools::tx_add_misincCount(DT) %>% txtools::tx_add_nucTotal()
    tmp <- round(DT$misincCount / DT$nucTotal, 6)
    tmp[DT$nucTotal < minReads] <- NA
    tibble::add_column(DT, diffToRef_Ratio = tmp)
}


# Start rate
add_StartRate <- function(DT, minCov = 50){
    tmp <- (DT$start_5p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, start_Ratio = tmp)
}

# Start rate
add_EndRate <- function(DT, minCov = 50){
    tmp <- (DT$end_3p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, end_Ratio = tmp)
}

# Start rate
add_StartRate_1bpDS <- function(DT, minCov = 50){
    tmp <- (DT$start_5p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, startRate_1bpDS = tmp) %>% txtools::tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$startRate_1bpDS <- c(utils::tail(DT$startRate_1bpDS, -1), NA)
        DT
    }) %>% txtools::tx_merge_DT()
    return(DT)
}

# 2Omescore calculation: Score A
detect_2OmeScore <- function(DT_low, DT_high, neighFlankSize = 5, minCov = 50){
    DT_low <- DT_low %>% add_StartRate(minCov = minCov)
    DT_high <- DT_high %>% add_StartRate(minCov = minCov)
    if(all(union(DT_low$gene, DT_high$gene) %in% intersect(DT_low$gene, DT_high$gene))){
        iGenes <- intersect(DT_low$gene, DT_high$gene)
    }else{stop("DT_low and DT_high do not share the same genes")}
    DTL1 <- DT_low %>% txtools::tx_split_DT()
    DTL2 <- DT_high %>% txtools::tx_split_DT()
    OUT <- lapply(iGenes, function(iG){
        RES <- data.table::data.table(DTL1[[iG]][,colnames(DTL1[[iG]]) %in% storm_baseCols, with = FALSE])
        tmpSR <- c(utils::tail(DTL1[[iG]]$start_Ratio / DTL2[[iG]]$start_Ratio, - 1), NA)
        sc_range <- (1 + neighFlankSize):(nrow(DTL1[[iG]]) - neighFlankSize)
        score <- rep(NA, nrow(DTL1[[iG]]))
        for(i in sc_range){
            score[i] <- log2(tmpSR[i]) - log2(mean(c(tmpSR[(i - neighFlankSize):(i - 1)],
                                                     tmpSR[(i + 1):(i + neighFlankSize)])))
        }
        RES$twoOme_score <- score
        return(RES)
    })
    do.call(what = rbind, OUT)
}

# Score A Birkedal et al., 2015
scoreA <- function(CHR, d = 6, minReads = 30){
    scoreA <- rep(NA, length(CHR)); names(scoreA) <- names(CHR)
    for(P in (d+1):(length(CHR)-d)){
        ePos <- CHR[P]
        leftFlank <- CHR[(P-d):(P-1)]
        rightFlank <- CHR[(P+1):(P+d)]
        if(sum(c(rightFlank, ePos, leftFlank)) < minReads){next()}
        m1 <- mean(leftFlank)
        m2 <- mean(rightFlank)
        sd1 <- stats::sd(leftFlank)
        sd2 <- stats::sd(rightFlank)
        scoreA[P] <- 1 - (((2 * ePos) + 1) /
                              ((0.5*abs(m1-sd1)) + ePos + (0.5*abs(m2-sd2)) + 1))
    }
    scoreA[scoreA < 0] <- 0
    return(scoreA)
}

# Adding score A to multiGene DT
add_scoreA_3p <- function(DT, minReads = 30, nCores = 1, flankSize = 6){
    DTL <- txtools::tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(x){
        tmp <- scoreA(x$end_3p, d = flankSize, minReads = minReads)
        tibble::add_column(x, scoreA_3p = tmp)
    }) %>% txtools::tx_merge_DT()
}

# Adding score A 5p_ends to multiGene DT
add_scoreA_5p <- function(DT, minReads = 30, nCores = 1, flankSize = 6){
    DTL <- txtools::tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(x){
        tmp <- scoreA(x$start_5p, minReads = minReads, d = flankSize)
        tmp <- c(utils::tail(tmp, -1), NA)
        tibble::add_column(x, scoreA_5p = tmp)
    }) %>% txtools::tx_merge_DT()
}

# Adding score C 3p_ends to multiGene DT
add_scoreC_3p <- function(DT, minReads = 30, nCores = 1, flankSize = 6){
    DTL <- txtools::tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(DT){
        tmp <- scoreC_3p(DT$end_3p, d = flankSize, minReads = minReads)
        tibble::add_column(DT, scoreC_3p = tmp)
    }) %>% txtools::tx_merge_DT()
}
# Start rate
add_StartRate <- function(DT, minCov = 50){
    tmp <- (DT$start_5p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, start_Ratio = tmp)
}

# Calculate deletion rate
deletion_rate <- function(DT, minReads = 30){
    DT <- txtools::tx_add_nucTotal(DT)
    tmp <- round(DT$`-` / DT$nucTotal, 6)
    tmp[DT$nucTotal < minReads] <- NA
    tibble::add_column(DT, deletion_Ratio = tmp)
}

# Calculates end rate 1bp downstream
add_EndRate1bpDS <- function(DT, minCov = 50){
    tmp <- (DT$end_3p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, endRate_1bpDS = tmp) %>% txtools::tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$endRate_1bpDS <- c(utils::tail(DT$endRate_1bpDS, -1), NA)
        DT
    }) %>% txtools::tx_merge_DT()
    return(DT)
}

# Calculates scoreZ
scoreZ <- function(CHR, d = 6, minReads = 30){
    score <- rep(NA, length(CHR)); names(score) <- names(CHR)
    for(P in (d+1):(length(CHR)-d)){
        ePos <- CHR[P]
        leftFlank <- CHR[(P-d):(P-1)]
        rightFlank <- CHR[(P+1):(P+d)]
        allpos = c(leftFlank, ePos, rightFlank)
        if(sum(allpos) > minReads){
            score[P] <- (ePos - mean(allpos)) / sd(allpos)
        }else{
            next()
        }
    }
    return(score)
}

# Calculates scoreZ on 3p
add_scoreZ_3p_util <- function(DT, minReads = 30, nCores = 1, flankSize = 6){
    DTL <- txtools::tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(x){
        tmp <- scoreZ(x$end_3p, d = flankSize, minReads = minReads)
        tibble::add_column(x, scoreZ_3p = tmp)
    }) %>% txtools::tx_merge_DT()
}

# Calculate scoreC on 3p
scoreC_3p <- function(CHR, d = 6, minReads = 30) {
    scoreC <- rep(NA, length(CHR)); names(scoreC) <- names(CHR)
    for(P in (d+1):(length(CHR)-d)){
        evalPos <- CHR[P]
        leftFlank <- CHR[(P - d):(P - 1)]
        rightFlank <- CHR[(P + 1):(P + d)]
        w = seq(1 - 0.1 * (d - 1), 1, 0.1)
        leftFlankNorm = leftFlank * w
        rightFlankNorm = rightFlank * rev(w)
        nreads=sum(leftFlank, evalPos, rightFlank)
        if (nreads > minReads) {
            scoreC[P] = 1 -
                (evalPos / (0.5*((sum(leftFlankNorm))/(sum(w)) +
                                     (sum(rightFlankNorm)) / (sum(rev(w))))))
        }else{
            next()
        }
    }
    return(scoreC)
}

# Scaling twoOme_score by nucleotide groups
scale_by_nuc <- function(DT){
    iNucs <- DT$refSeq %>% unique()
    for(iN in iNucs){
        selPos <- which(DT$refSeq == iN)
        DT$twoOme_score[selPos] <- scale(DT$twoOme_score[selPos])
    }
    return(DT)
}
