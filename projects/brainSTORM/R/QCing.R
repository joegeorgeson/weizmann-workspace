plot_metricsNArate <- function(STORM){
    sets <- names(STORM$CALLS)
    tmpL <- lapply(sets, function(set){
        RNAmods <- names(STORM$CALLS[[set]])
        lapply(RNAmods, function(RNAmod_i){
            tmpCALLS <- STORM$CALLS[[set]][[RNAmod_i]]
            tmpC_met <- tmpCALLS[,intersect(metricsList[[RNAmod_i]], colnames(tmpCALLS)), with = FALSE]
            tmpO <- tmpC_met[tmpCALLS$refSeq %in% RNAMod_baseNuc[[RNAmod_i]],] %>% is.na %>% colMeans() %>% round(3)
            data.table(set = set, RNAmod = RNAmod_i, metric = names(tmpO), NA_rate = tmpO)
        }) %>% do.call(what = "rbind")
    }) %>% do.call(what = "rbind")
    ggplot(tmpL, aes(x = metric, y = NA_rate)) +
        geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        facet_grid(set~RNAmod) + coord_flip()
}

plot_predsNApct <- function(STORM){
    tmpO <- qc_predsNApct(STORM)
    ggplot(tmpO, aes(x = RNAmod, y = NApred_pct)) +
        geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        facet_grid(set~.) + coord_flip() + ylab("Prediction NA %")
}

#' QC Predictions NA percentage
#'
#' @param STORM
#'
#' @return
#' @export
#'
qc_predsNApct <- function(STORM){
    if("CALLS" %in% names(STORM)){
        sets <- names(STORM$CALLS)
        tmpO <- lapply(sets, function(set_i){
            RNAmods <- names(STORM$CALLS[[set_i]])
            lapply(RNAmods, function(RNAmod_i){
                tmpDT <- STORM$CALLS[[set_i]][[RNAmod_i]]
                naRate <- tmpDT[tmpDT$refSeq %in% RNAMod_baseNuc[[RNAmod_i]],]$pred %>%
                    is.na() %>% mean() %>% multiply_by(100) %>% round(2)
                data.table(set = set_i, RNAmod = RNAmod_i, NApred_pct = naRate)
            }) %>% do.call(what = "rbind")
        }) %>% do.call(what = "rbind")
        tmpO
    }else{stop("CALLS is not an element of STORM object. Use ",
               "storm_assignScores() and call_RNAmods_logRegMods() to compute ",
               "RNAmod predictions")}
}
