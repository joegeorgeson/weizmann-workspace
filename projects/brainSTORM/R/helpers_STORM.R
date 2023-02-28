# Magrittr Pipe Operator
`%>%` <- magrittr::`%>%`

# Make Temporary directory for STORMseq processing
mkTmpDir <- function(){
    if(!dir.exists("./STORMtmp_dir")){dir.create("./STORMtmp_dir")}
}

#' Remove STORM temporary directory
#'
#' Removes temporary directory "./STORMtmp_dir" created for alignment steps.
#'
#' @return
#' @export
#'
#' @examples
rmTmpDir <- function(){
    if(dir.exists("./STORMtmp_dir")){unlink("./STORMtmp_dir", recursive = TRUE)}
}

# Extract gene sequences from genome and geneAnnotation
getGeneSeqsfromGenome <- function(geneAnnot, genome, nCores = 1){
    txtools:::check_mc_windows(nCores)
    txtools:::check_GA_genome_chrCompat(geneAnnot = geneAnnot, genome = genome)
    parallel::mclapply(mc.cores = nCores, seq_along(geneAnnot), function(i){
        selGene <- geneAnnot[i]
        iGene <- as.character(selGene$name)
        iChr <- as.character(GenomicRanges::seqnames(selGene))
        iStr <- as.character(selGene@strand)
        iGA <- selGene
        iBlocks <- S4Vectors::mcols(iGA)$blocks %>% txtools:::if_IRangesList_Unlist() %>%
            IRanges::shift(IRanges::start(iGA) - 1)
        SEQ <- stringr::str_sub(genome[[iChr]], start = IRanges::start(iBlocks),
                                end = IRanges::end(iBlocks)) %>% paste(collapse = "") %>%
            Biostrings::DNAString()
        if(iStr == "-") {
            SEQ <- Biostrings::reverseComplement(SEQ)
        }
        SEQ
    }) %>% Biostrings::DNAStringSet()
}

# Add results to a STORM object. Remove scores if metric is already present.
hlpr_add_REScols <- function(STORM_RES, REScols){
    iMetric <- unique(REScols[,"metric"]) %>% as.character()
    # remove results for identical metric
    if("metric" %in% names(STORM_RES)){
        STORM_RES <- STORM_RES[STORM_RES$metric != iMetric,]
    }
    STORM_RES <- rbind(STORM_RES, REScols)
    STORM_RES
}

#' Check RNAmods
#'
#' Check which RNAmods are possible to predict based on RNAmod metrics present
#' in STORM object
#'
#' @param STORM
#'
#' @return character
#' @export
check_whichRNAmods <- function(STORM){
    if("RES" %in% names(STORM)){
        tmp <- lapply(brainSTORM:::metricsList, function(x) any(unique(STORM$RES$metric) %in% x)) %>%
            unlist()
        return(names(tmp[tmp]))
    }else{stop("STORM object does not include 'RES' element")}
}

# Confusion matrix helper
hlpr_confusMatOutcome <- function(x){
    x$outcome <- paste0(ifelse(x$truth == x$pred, "T", "F"),
                        ifelse(as.logical(x$pred), "P", "N"))
    x$outcome[is.na(x$truth) | is.na(x$pred)] <- "NA"
    y <- x[x$outcome == "NA",]
    NAsum <- tapply(y$Freq, y$refSeq, "sum")
    out <- rbind(x[x$outcome %in% c("TN", "TP", "FP", "FN") & !is.na(x$refSeq),],
          data.table::data.table(truth = NA, pred = NA,
                                 refSeq = names(NAsum),
                                 Freq = as.numeric(NAsum),
                                 outcome = "NA"))
    return(out)
}

hlp_removeColumnIfPresent <- function(DT, colName){
    if(colName %in% colnames(DT)){
        DT[, colnames(DT) != colName, with = FALSE]
    }else{DT}
}

check_DT <- function(DT){
    if(!data.table::is.data.table(DT)){
        if(!is.data.frame(DT)){
            stop("DT must be either a data.table or a data.frame")
        }else{
            if(is.data.frame(DT)){DT <- data.table::data.table(DT)}
        }
    }
    return(DT)
}
