storm_summary <- function(STORM){
    RNAmods <- names(STORM$CALLS[[1]])
    OUT <- lapply(RNAmods, function(RNAmod_i){
        tmp_out <- STORM$DATA[[1]][,storm_baseCoorCols, with = FALSE]
        SETS <- levels(STORM$META$set)
        if(is.null(STORM$CALLS[[SETS[1]]][[RNAmod_i]]$logist_Score)){
            warning("logist_Score was not found for ", RNAmod_i, ". Results will",
                    " not be summarized in STORM$SUMMARY. Check if metrics where ",
                    "placed in STORM$CALLS[[1]]$", RNAmod_i)
            return(NULL)
        }else{
            tmp_log <- lapply(SETS, function(x){
                STORM$CALLS[[x]][[RNAmod_i]]$logist_Score
            }) %>% do.call(what = cbind) %>% magrittr::set_colnames(paste("logScore",
                                                                          SETS, sep = "_"))
            tmp_lin <- lapply(SETS, function(x){
                STORM$CALLS[[x]][[RNAmod_i]]$linear_Score
            }) %>% do.call(what = cbind) %>% magrittr::set_colnames(paste("linScore",
                                                                          SETS, sep = "_"))
            tmp_prd <- lapply(SETS, function(x){
                STORM$CALLS[[x]][[RNAmod_i]]$pred
            }) %>% do.call(what = cbind) %>% magrittr::set_colnames(paste("pred",
                                                                          SETS, sep = "_"))
            tmp_out <- cbind(tmp_out, tmp_lin, tmp_log, tmp_prd)
            return(tmp_out)
        }
    }) %>% magrittr::set_names(RNAmods)
    STORM$SUMMARY <- OUT
    STORM
}


# Tables of files from FASTQ to expected BAM and RDS targets
files_table <- function(META){
    fastq <- META$FASTQ
    bam <- META$BAM
    lce <- gsub(pattern = "Aligned.out.sorted.bam",
                replacement = "Aligned.out.sorted.lce.txt", META$BAM)
    rds <- gsub(pattern = "Aligned.out.sorted.bam",
                replacement = "Aligned.out.sorted.txDT.rds", META$BAM)
    # Table
    tmpDT <- data.table::data.table(FASTQ = c(fastq),
                                    BAM = c(bam),
                                    BAM_ok = file.exists(c(bam)),
                                    lce = c(lce),
                                    lce_ok = file.exists(c(lce)),
                                    rds = c(rds),
                                    rds_ok = file.exists(c(rds)))
    return(tmpDT)
}

# FASTQ Duplication rate (library complexity)
fastq_dupRate <- function(FASTQs_pahts, nCores){
    parallel::mclapply(mc.cores = nCores, FASTQs_pahts, function(file){
        tmp <- ShortRead::readFastq(file)
        tmp2 <- ShortRead::readFastq(gsub(file, pattern = "R1", replacement = "R2"))
        dupR <- paste(tmp@sread, tmp2@sread, sep = "") %>% duplicated() %>% mean
        return(dupR)
    }) %>% unlist
}

#' Nucleotide frequency report
#'
#' @param META
#' @param nCores
#' @param firstN
#' @param pairedEnd
#'
#' @return
#' @export
#'
#' @examples
fastq_nucFreq <- function(META, nCores, firstN = 1e4, pairedEnd = TRUE){
    parallel::mclapply(mc.cores = nCores, META$FASTQ, function(file){
        tmp <- readLines(file, firstN * 4)[seq(2, firstN*4, 4)] %>% Biostrings::DNAStringSet()
        r2File <- gsub(file, pattern = "R1", replacement = "R2")
        if(pairedEnd){
            tmp2 <- readLines(r2File, firstN * 4)[seq(2, firstN*4, 4)] %>%
                Biostrings::DNAStringSet() %>% Biostrings::complement()
            mR1R2 <- paste(tmp, tmp2, sep = "")
            nucFreq <- mR1R2 %>% stringr::str_split(pattern = "") %>% unlist %>% table
        }else if(!pairedEnd){
            nucFreq <- tmp %>% stringr::str_split(pattern = "") %>% unlist %>% table
        }
        return(nucFreq[c("A", "C", "G", "T")])
    }) %>% do.call(what = cbind) %>% magrittr::set_colnames(META$id)
}

#' nucleotide frequency plot
#'
#' @param nucF_x
#' @param subtitle
#'
#' @return
#' @export
#'
#' @examples
gg_nucFreq <- function(nucF_x, subtitle){
    tmp <- prop.table(nucF_x, margin = 2) %>% data.frame %>%
        tibble::rownames_to_column(var = "nuc") %>%
        tidyr::pivot_longer(cols = -nuc, names_to = "Sample", values_to = "Ratio")
    ggplot2::ggplot(tmp, ggplot2::aes(x = Sample, y = Ratio, fill = nuc)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_brewer(palette="Set1") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle("Nucleotide Frequency per library", subtitle = subtitle)
}

#' Alignment report table
#'
#' @param META
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
reads_report <- function(META, nCores = 4){
    DT <- files_table(META)
    if(all(DT$BAM_ok) & all(DT$rds_ok)){
        res1 <- parallel::mclapply(mc.cores = nCores, DT$FASTQ, function(x){
            tmp <- ShortRead::readFastq(x)
            length(tmp)
        }) %>% unlist
        res2 <- parallel::mclapply(mc.cores = nCores, DT$BAM[DT$BAM_ok], function(x){
            tmp2 <- Rsamtools::scanBam(x)
            tmp2[[1]]$qname %>% unique %>% length
        }) %>% unlist
        res3 <- parallel::mclapply(mc.cores = nCores, DT$rds[DT$rds_ok], function(x){
            tmplog <- data.table::fread(gsub(pattern = "rds", "log", x), header = F)
            tmplog[grep(tmplog$V1, pattern =  "unique reads"),]$V2 %>% as.numeric()
        }) %>% unlist
        tmpDT <- data.table::data.table(sample = META$id,
                                        FASTQ_reads = res1,
                                        BAM_aligns = NA,
                                        tx_starts = NA)
        tmpDT$BAM_aligns[DT$BAM_ok] <- res2
        tmpDT$tx_starts[DT$rds_ok] <- res3
        tmpDT$pC_BAM <- round(tmpDT$BAM_aligns / tmpDT$FASTQ_reads * 100, 2)
        tmpDT$pC_tx <- round(tmpDT$tx_starts / tmpDT$BAM_aligns * 100, 2)
        return(tmpDT)
    }else{
        stop("Some files are missing.")
    }
}

#' Alignment report plots
#'
#' @param rReport
#' @param species
#'
#' @return
#' @export
#'
#' @examples
gg_readStats <- function(rReport, species){
    tmpDT <- rReport[, -c(5, 6)]
    tmpDT$FASTQ <- tmpDT$FASTQ_reads - tmpDT$BAM_aligns
    tmpDT$BAM <- tmpDT$BAM_aligns - tmpDT$tx_starts
    tmpDT$txDT <- tmpDT$tx_starts
    tmpDT$BAM[is.na(tmpDT$BAM)] <- 0; tmpDT$txDT[is.na(tmpDT$txDT)] <- 0
    tmpDT <- tidyr::pivot_longer(data = tmpDT[,c("sample", "FASTQ", "BAM", "txDT")],
                                 cols = c("FASTQ", "BAM", "txDT"), names_to = "Reads")
    tmpDT$Reads <- factor(tmpDT$Reads, levels = c("FASTQ", "BAM", "txDT"))
    tmpDT$value <- magrittr::divide_by(tmpDT$value, 1e6)
    t_GG1 <- ggplot2::ggplot(tmpDT, ggplot2::aes(x = sample, y = value, fill = Reads)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(paste("Reads in", species, "libraries by processing step"),) +
        ggplot2::ylab("Million reads") +
        ggplot2::xlab("Samples")
    # Proportion
    tmpDT <- rReport[, -c(5, 6)]
    tmpDT$FASTQ <- tmpDT$FASTQ_reads - tmpDT$BAM_aligns
    tmpDT$BAM <- tmpDT$BAM_aligns - tmpDT$tx_starts
    tmpDT$txDT <- tmpDT$tx_starts
    tmpDT$BAM[is.na(tmpDT$BAM)] <- 0; tmpDT$txDT[is.na(tmpDT$txDT)] <- 0
    tmpDT <- data.frame(sample = tmpDT$sample,
                        apply(tmpDT[,c("FASTQ", "BAM", "txDT")], 1, prop.table) %>% t)
    tmpDT <- tidyr::pivot_longer(data = tmpDT[,c("sample", "FASTQ", "BAM", "txDT")],
                                 cols = c("FASTQ", "BAM", "txDT"), names_to = "Reads")
    tmpDT$Reads <- factor(tmpDT$Reads, levels = c("FASTQ", "BAM", "txDT"))

    t_GG2 <- ggplot2::ggplot(tmpDT, ggplot2::aes(x = sample, y = value, fill = Reads)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(paste("Proportion of reads in", species, "libraries by processing step"),) +
        ggplot2::ylab("Proportion") +
        ggplot2::xlab("Samples")
    return(list(t_GG1, t_GG2))
}

# Alignment and transcript data processing efficiency
ggAlignEffPlot <- function(META, rReport){
    tmp <- cbind(META, rReport[,-1])
    tmpGG1 <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$libTreat, y = tmp$pC_BAM, colour = tmp$libTreat)) +
        ggplot2::geom_boxplot() + ggplot2::geom_point(colour = "black") +
        ggplot2::ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ggplot2::ylab("% Aligned reads") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    tmpGG2 <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$bioTreat, y = tmp$pC_BAM, colour = tmp$bioTreat)) +
        ggplot2::geom_boxplot() + ggplot2::geom_point(colour = "black") +
        ggplot2::ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ggplot2::ylab("% Aligned reads") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    tmpGG3 <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$RTase, y = tmp$pC_BAM, colour = tmp$RTase)) +
        ggplot2::geom_boxplot() + ggplot2::geom_point(colour = "black") +
        ggplot2::ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ggplot2::ylab("% Aligned reads") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    list(tmpGG1, tmpGG2, tmpGG3)
}

#' Library complexity extrapolation plot
#'
#' @param META
#' @param tab_name
#' @param speciesName
#'
#' @return
#' @export
#'
#' @examples
gg_lce <- function(META, tab_name, speciesName = ""){
    lceFiles <- gsub(META$BAM, pattern = ".bam$", replacement = ".lce.txt") %>%
        magrittr::set_names(META$id)
    if(!all(file.exists(lceFiles))){
        stop("Report files missing:\n", paste(lceFiles[!file.exists(lceFiles)], collapse = " \n"))
    }
    tmp <- lapply(lceFiles, function(x) data.table::fread(x)) %>%
        magrittr::set_names(META$id)
    tmp <- lapply(names(tmp), function(x){
        cbind(tmp[[x]], id = x)
    }) %>% do.call(what = rbind) %>% data.table::data.table()
    data.table::fwrite(x = tmp, file = tab_name, sep = "\t")
    cat("Library complexity table output:", tab_name)
    tmpT <- table(tmp$TOTAL_READS)
    e_reads <- names(tmpT)[tmpT == max(tmpT)] %>% as.numeric %>% max
    tmp <- tmp[tmp$TOTAL_READS == e_reads,]
    ggOUT <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$id, y = tmp$EXPECTED_DISTINCT)) +
        ggplot2::geom_bar(stat="identity", color="black", position= ggplot2::position_dodge()) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin= tmp$LOWER_0.95CI, ymax= tmp$UPPER_0.95CI), width=.2,
                               position= ggplot2::position_dodge(1)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(label = paste("Library complexity extrapolation -", speciesName),
                         sub = paste("Expected distinct reads at", e_reads, "depth. CI = 95%")) +
        ggplot2::ylab("Expected distinct reads") +
        ggplot2::xlab("Samples")
    return(ggOUT)
}

# is.Modified functions ##########################################################

# RNAmod logical vectors from character or factor
is.pseudoU <- function(nucs){nucs %in% c("Y", "Ym")}
is.2Ome <- function(nucs){nucs %in% c("Am", "Gm", "Um", "Cm", "Ym")}
is.m5C <- function(nucs){nucs %in% c("m5C")}
is.ac4C <- function(nucs){nucs %in% c("ac4C")}
is.m1A <- function(nucs){nucs %in% c("m1A")}
is.m7G <- function(nucs){nucs %in% c("m7G")}
is.m3U <- function(nucs){nucs %in% c("m3U")}
is.Am <- function(nucs){nucs %in% c("Am")}
is.Cm <- function(nucs){nucs %in% c("Cm")}
is.Gm <- function(nucs){nucs %in% c("Gm")}
is.Tm <- function(nucs){nucs %in% c("Um", "Ym")}
is.m62A <- function(nucs){nucs %in% c("m62A")}
is.m1acp3Y <- function(nucs){nucs %in% c("m1acp3Y")}
is.m6A <- function(nucs){nucs %in% c("m6A")}

# Modeling functions #####

#  Balanced groups assignment of two level vectors for K-fold cross validation
CV_balancedGroups <- function(x, k){
    x <- factor(x)
    blocks <- seq(1, k)
    groups <- rep(NA, length(x))
    tmpLog1 <- x == levels(x)[1]
    tmpLog2 <- x == levels(x)[2]
    groups[tmpLog1] <- sample(rep(blocks, ceiling(sum(tmpLog1)/length(blocks))))[seq(groups[tmpLog1])]
    groups[tmpLog2] <- sample(rep(blocks, ceiling(sum(tmpLog2)/length(blocks))))[seq(groups[tmpLog2])]
    groups
}
# TODO: Make a generalized K-fold balanced groups sampler

# Extract metrics from k-CFV res
extract_AUC_sen_spe <- function(res, RNAmod, strategy, organism){
    data.frame(organism = organism,
               strategy = strategy,
               RNAmod = RNAmod,
               AUC = unlist(lapply(seq(res), function(i) res[[i]][RNAmod, "AUC"])),
               sensitivity = unlist(lapply(seq(res), function(i) res[[i]][RNAmod, "Sensitivity"])),
               specificity = unlist(lapply(seq(res), function(i) res[[i]][RNAmod, "Specificity"])))
}

# Matthew Correlation Coefficient
matthewCorrCoeff <- function(scores, successes){
    MCC <- NULL
    known <- rep(0, length(scores)); known[successes] <- 1
    NAs <- which(is.na(scores))
    scores[is.na(scores)] <- 0
    for (i in 1:length(scores)){
        thr <- scores[i]
        calls <- as.numeric(scores >= thr)
        TP <-  sum(calls & known)
        FP <-  sum(calls) - TP
        TN <- sum(!calls & !known)
        FN <- sum(!calls) - TN
        MCC[i] <- (TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    }
    MCC[is.nan(MCC)] <- 0
    MCC[is.infinite(MCC)] <- 0
    MCC[NAs] <- NA
    return(MCC)
}

# MCC, using threshold removing NA scores
MCC <- function(scores, successes, thr){
    successes <- successes[!is.na(scores)]
    scores <- scores[!is.na(scores)]
    calls <- as.numeric(scores >= thr)
    known <- rep(0, length(scores)); known[successes] <- 1
    TP <-  sum(calls & known)
    FP <-  sum(calls) - TP
    TN <- sum(!calls & !known)
    FN <- sum(!calls) - TN
    MCC <- (TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if(is.nan(MCC)){MCC <- 0}
    MCC
}

# All combinations of strings in character vector
all_comb_allNelements <- function(x){
    lapply(seq_along(x), function(y) utils::combn(x, y, simplify = FALSE)) %>%
        unlist(recursive = FALSE)
}

#' Logistic Model training using STORM$CALLS object.
#'
#' @param STORM
#' @param isModFun
#' @param RNAmod
#' @param accNArate
#' @param thresholds
#' @param minMCC
#'
#' @return
#' @export
trainLogModel_RNAmods <- function(STORM, isModFun, RNAmod, thresholds,
                                  accNArate = TRUE, minMCC = 0.7){
    CALLS <- lapply(STORM$CALLS, function(x) data.frame(x[[RNAmod]])) %>% do.call(what = "rbind")
    CALLS$isMod <- isModFun(CALLS$nuc)
    varNames <- names(CALLS)[names(CALLS) %in% allMetrics]
    allCombVar <- all_comb_allNelements(varNames)
    combNames <- paste0("Comb_", seq_along(allCombVar))
    OUT <- lapply(allCombVar, function(selVars){
        tData <- CALLS[CALLS$refSeq %in% RNAMod_baseNuc[[RNAmod]], c(selVars, "isMod")]
        tData$isMod <- as.numeric(tData$isMod)
        logMod <- stats::glm(isMod ~ ., family = stats::binomial(), data = tData) %>% suppressWarnings()
        CALLS$logPred <- stats::predict(logMod, newdata = CALLS, type = "response")
        CALLS <- CALLS[CALLS$refSeq %in% RNAMod_baseNuc[[RNAmod]],]
        return(list(MCCs = sapply(thresholds, function(x) MCC(CALLS$logPred, CALLS$isMod, x)),
                    NArate = mean(is.na(CALLS$logPred))))
    })
    MCCmat <- lapply(OUT, function(x) x$MCCs) %>% do.call(what = cbind) %>%
        magrittr::set_rownames(thresholds) %>% magrittr::set_colnames(combNames)
    NArate <- lapply(OUT, function(x) x$NArate) %>% unlist()
    medianMCC <- apply(MCCmat, 2, "median")
    maxMCC <-  apply(MCCmat, 2, "max")
    if(accNArate){
        medianMCC <- medianMCC
        maxMCC <- maxMCC
        modScore <- ((1 - NArate) * medianMCC)
        modScore2 <- ((1 - NArate) * maxMCC)
        selComb <- which.max(modScore * modScore2)
        selThr <- which.max(MCCmat[,selComb])
    }else{
        selComb <- which.max(medianMCC)
        selThr <- which.max(MCCmat[,selComb])
    }
    if(MCCmat[selThr, selComb] < minMCC){
        warning("Model for RNAmod:", RNAmod, " was set to NULL as obtained MCC:",
                round(MCCmat[selThr, selComb], 4), " is less than minimum MCC:",
                minMCC)
        return(list(logiMod = NULL, thr = thresholds[selThr], vars = NULL,
                    MCCmat = MCCmat, MCC = NULL, NArate = NULL,
                    onNucs = RNAMod_baseNuc[[RNAmod]], medianMCC = medianMCC,
                    maxMCC = maxMCC, NArates = NArate, selComb = NULL))
    }
    # FinalModel
    tData <- CALLS[CALLS$refSeq %in% RNAMod_baseNuc[[RNAmod]], c(allCombVar[[selComb]], "isMod")]
    tData$isMod <- as.numeric(tData$isMod)
    logMod <- stats::glm(isMod ~ ., family = stats::binomial(), data = tData) %>% suppressWarnings()
    return(list(logiMod = logMod, thr = thresholds[selThr], vars = allCombVar[[selComb]],
                MCCmat = MCCmat, MCC = MCCmat[selThr, selComb], NArate = NArate[selComb],
                onNucs = RNAMod_baseNuc[[RNAmod]], medianMCC = medianMCC, maxMCC = maxMCC,
                NArates = NArate, selComb = selComb))
}

#' Train all logistic models
#'
#' @param STORM list. STORM object
#' @param thrRange numeric vector. set of thresholds used to search for a maximum
#' MatthewsCorrelationCoefficient
#' @param accNArate logical. Account for NA rate for model selection
#' @param minMCC numeric. Minimum MCC needed for RNAmod model to be output.
#' Set to -1 to force all models to be output.
#'
#' @return
#' @export
trainLogModel_allRNAmods <- function(STORM, thrRange = seq(-0.1, 0.5, 0.01),
                                     accNArate = TRUE, minMCC = 0.7){
    RNAmods <- check_whichRNAmods(STORM)
    lapply(RNAmods, function(RNAmod_i){
        trainLogModel_RNAmods(STORM = STORM, isModFun = RNAModFunList[[RNAmod_i]],
                              RNAmod = RNAmod_i, thresholds = thrRange,
                              accNArate = accNArate, minMCC = minMCC)
    }) %>% magrittr::set_names(RNAmods)
}


#' Logistic regressions scores
#'
#' Update STORM$CALLS with logistic and linear scores and prediction
#'
#' @param STORM
#' @param logMods
#'
#' @return
#' @export
call_RNAmods_logRegMods <- function(STORM, logMods, nucBias = TRUE){
    STORM$CALLS <- lapply(names(STORM$CALLS), function(set){
        RNAmods <- intersect(names(STORM$CALLS[[set]]), names(logMods))
        tmp <- lapply(RNAmods, function(RNAmod){
            if(is.null(logMods[[RNAmod]]$logiMod)){
                warning("There is no model for RNAmod:", RNAmod,
                        " predictions and logistic",
                        " scores, were not calculated.")
                return(STORM$CALLS[[set]][[RNAmod]])
            }
            if(all(logMods[[RNAmod]]$vars %in% colnames(STORM$CALLS[[set]][[RNAmod]]))){
                nData <- data.frame(STORM$CALLS[[set]][[RNAmod]])
                STORM$CALLS[[set]][[RNAmod]]$logist_Score <-
                    stats::predict.glm(logMods[[RNAmod]]$logiMod, newdata = nData, type = "response")
                STORM$CALLS[[set]][[RNAmod]]$linear_Score <-
                    stats::predict.glm(logMods[[RNAmod]]$logiMod, newdata = nData, type = "link")
                STORM$CALLS[[set]][[RNAmod]]$pred <-
                    STORM$CALLS[[set]][[RNAmod]]$logist_Score >= logMods[[RNAmod]]$thr
                if(nucBias){
                    notBaseNuc <- !(STORM$CALLS[[set]][[RNAmod]]$refSeq %in% logMods[[RNAmod]]$onNucs)
                    STORM$CALLS[[set]][[RNAmod]]$logist_Score[notBaseNuc] <- NA
                    STORM$CALLS[[set]][[RNAmod]]$linear_Score[notBaseNuc] <- NA
                    STORM$CALLS[[set]][[RNAmod]]$pred[notBaseNuc] <- NA
                }
                STORM$CALLS[[set]][[RNAmod]]$pred[is.na(STORM$CALLS[[set]][[RNAmod]]$logist_Score)] <- NA
                return(STORM$CALLS[[set]][[RNAmod]])
            }else{
                warning("Some model metrics are not found in STORM$CALLS ",
                        "for RNAmod:", RNAmod, " predictions and logistic",
                        " scores, were not calculated.")
                return(STORM$CALLS[[set]][[RNAmod]])
            }
            return(tmp)
        }) %>% magrittr::set_names(RNAmods)
    }) %>% magrittr::set_names(names(STORM$CALLS))
    STORM
}

# Save as RDS if object is not already RDS
saveAsRDSIfNotAlready <- function(object, fileName){
    if(!file.exists(fileName)){saveRDS(object, fileName)}
}

# Has duplicates?
hasDups <- function(x){
    sum(duplicated(x)) > 0
}

# Limit STORM object to genes in geneAnnot argument
storm_reduceToGA <- function(STORM, geneAnnot){
    STORM$DATA <- lapply(STORM$DATA, function(DT) DT[gene %in% geneAnnot$name,])
    return(STORM)
}

# storm_calls: Makes predictions based on RF and CutPointR for each modification
storm_calls <- function(STORM, RNAmodList, RF_list, CP_models){
    RNAmods_vec <- names(RNAmodList)
    coorSys <- STORM$RES[,names(STORM$RES) %in% c("chr", "gencoor", "strand",
                                                  "gene", "txcoor", "pos",
                                                  "refSeq", "nuc"), with = FALSE]
    iSets <- levels(STORM$META$set)
    STORM$CALLS <- lapply(seq_along(RNAmods_vec), function(i){
        pred_RES <- lapply(seq_along(iSets), function(j){
            lapply(CP_models[[i]], function(CPmod){
                tmpPat <- paste0("^", gsub("_TGIRT", "", CPmod$predictor))
                selVar <- grepl(pattern = tmpPat, names(STORM$RES)) &
                    grepl(pattern = iSets[j], names(STORM$RES))
                tmpDat <- data.frame(x = STORM$RES[[which(selVar)]])
                CPmod$predictor <- "x"
                as.logical(stats::predict(CPmod, tmpDat))
            }) %>% do.call(what = data.frame) %>% unname %>% apply(1, all)
        }) %>% do.call(what = cbind) %>% data.table::data.table() %>%
            magrittr::set_names(paste("pred", RNAmods_vec[i], iSets, sep = "_"))

        scor_RES <- lapply(seq_along(iSets), function(j){
            lapply(CP_models[[i]], function(CPmod){
                tmpPat <- paste0("^", gsub("_TGIRT", "", CPmod$predictor))
                selVar <- grepl(pattern = tmpPat, names(STORM$RES)) &
                    grepl(pattern = iSets[j], names(STORM$RES))
                STORM$RES[,which(selVar), with = FALSE]
            }) %>% do.call(what = cbind)
        }) %>% do.call(what = cbind) %>% data.table::data.table()

        tree_RES <- lapply(seq_along(iSets), function(j){
            RF_mod <- RF_list[[i]]
            RF_vars <- names(RF_mod$forest$xlevels)
            newDat <- lapply(RF_vars, function(RF_v){
                tmpP <- gsub(pattern = "_TGIRT", replacement = "", x = RF_v)
                selVar <- grepl(pattern = tmpP, x = names(STORM$RES))
                selSmp <- grepl(pattern = iSets[j], x = names(STORM$RES)) & selVar
                newDat <- data.frame(x = STORM$RES[, selSmp, with = FALSE]) %>% magrittr::set_names(RF_v)
            }) %>% do.call(what = cbind)
            newDat[is.na(newDat)] <- 0
            out <- stats::predict(RF_mod, newDat, norm.votes = TRUE, type = "vote")
            unname(out[,"TRUE"])
        }) %>% do.call(what = cbind) %>% data.table::data.table() %>%
            magrittr::set_names(paste("votes", RNAmods_vec[i], iSets, sep = "_"))

        data.table::data.table(coorSys, pred_RES, tree_RES, scor_RES)
    }) %>% magrittr::set_names(RNAmods_vec)
    return(STORM)
}


storm_confMatrix <- function(STORM, nucBias = TRUE){
    if(!"CALLS" %in% names(STORM)){
        stop("CALLS is not an element of STORM, use storm_assignScores() to ",
             "assign metrics to the")
    }
    lapply(names(STORM$CALLS), function(set_i){
        lapply(names(STORM$CALLS[[set_i]]), function(RNAmod_i){
            RNAModFun <- RNAModFunList[[RNAmod_i]]
            selCALLS <- STORM$CALLS[[set_i]][[RNAmod_i]]
            if(nucBias){
                selCALLS <- selCALLS[selCALLS$refSeq %in% RNAMod_baseNuc[[RNAmod_i]],]
            }
            if(is.null(selCALLS$pred)){
                return(NULL)
            }
            tmpT <- table(truth = factor(RNAModFun(selCALLS$nuc), levels = c("TRUE", "FALSE")),
                          pred = factor(selCALLS$pred, levels = c("TRUE", "FALSE")),
                          refSeq = selCALLS$refSeq, useNA = "always") %>%
                data.frame() %>% data.table::data.table() %>% hlpr_confusMatOutcome()
            return(data.table(outcome = tmpT$outcome,
                              refSeq = tmpT$refSeq,
                              RNAmod = RNAmod_i,
                              set = set_i,
                              freq = tmpT$Freq))
        }) %>% do.call(what = rbind)
    }) %>% do.call(what = rbind)
}

storm_confMatrix4plot <- function(STORM){
    lapply(names(STORM$CALLS), function(set_i){
        lapply(names(STORM$CALLS[[set_i]]), function(RNAmod_i){
            RNAModFun <- RNAModFunList[[RNAmod_i]]
            tmpT <- table(truth = STORM$CALLS[[set_i]][[RNAmod_i]]$nuc %>% RNAModFun(),
                          pred = STORM$CALLS[[set_i]][[RNAmod_i]]$pred,
                          refSeq = STORM$CALLS[[set_i]][[RNAmod_i]]$refSeq) %>%
                data.frame() %>% data.table::data.table() %>% hlpr_confusMatOutcome()
            data.table(outcome = rep(tmpT$outcome, tmpT$Freq),
                       refSeq = rep(tmpT$refSeq, tmpT$Freq),
                       RNAmod = RNAmod_i,
                       set = set_i)
        }) %>% do.call(what = rbind)
    }) %>% do.call(what = rbind)
}

#' Evaluate predictions
#'
#' Evaluate predictions using known modified nucleotides identity
#'
#' @param STORM list. STORM object
#' @param nucBias logical. Nucleotide bias: only cases in the base nucleotide
#' of the respective RNAmod are considered.
#' @param roundDig numeric. Number of digits to round the resulting specificity,
#' sensitivity, FDR, and query rate values.
#'
#' @return
#' @export
#'
storm_evalPreds <- function(STORM, nucBias = TRUE, roundDig =4){
    confMat <- storm_confMatrix(STORM, nucBias = nucBias)
    lapply(unique(confMat$set), function(set_i){
        lapply(unique(confMat$RNAmod), function(RNAmod_i){
            RNAModFun <- RNAModFunList[[RNAmod_i]]
            tmpC <- confMat[confMat$RNAmod == RNAmod_i & confMat$set == set_i,]
            tmpC <- tapply(tmpC$freq, tmpC$outcome, "sum")
            sens <- tmpC["TP"] / (tmpC["TP"] + tmpC["FN"])
            spec <- tmpC["TN"] / (tmpC["FP"] + tmpC["TN"])
            FDR <- tmpC["FP"] / (tmpC["FP"] + tmpC["TP"])
            nMod <- sum(RNAModFun(STORM$CALLS[[set_i]][[RNAmod_i]]$nuc))
            nPred <- sum(tmpC["TP"], tmpC["TN"], tmpC["FP"], tmpC["FN"], na.rm =  TRUE)
            nTrial <- sum(tmpC["TP"], tmpC["TN"], tmpC["FP"], tmpC["FN"], tmpC["NA"], na.rm =  TRUE)
            data.frame(set = set_i,
                       RNAmod = RNAmod_i,
                       TP = tmpC["TP"],
                       FP = tmpC["FP"],
                       TN = tmpC["TN"],
                       FN = tmpC["FN"],
                       sens = round(sens, roundDig),
                       spec = round(spec, roundDig),
                       FDR = round(FDR, roundDig),
                       n_knownRNAmod = nMod,
                       queryRate = round(nPred / nTrial, roundDig)
            )
        }) %>% do.call(what = "rbind")
    }) %>% do.call(what = "rbind") %>% magrittr::set_rownames(NULL) %>%
        data.table::data.table()
}


#' Restore txDT Genomic Coordinate System
#'
#' @param txDT data.table
#' @param geneAnnot GenomicRanges Gene annotation loaded via tx_load_bed()
#' @param nCores integer
#'
#' @return data.table
#' @export
restoreGenCoors <- function (txDT, geneAnnot, nCores = 1){
    txtools:::check_mc_windows(nCores)
    if (all(txDT$gene %in% geneAnnot$name)) {
        txLengths <- txtools::tx_get_geneLengths(txDT)
        genCoorSys <- parallel::mclapply(mc.cores = nCores,
                                         as.character(unique(txDT$gene)),
                                         function(iGene){
                                             tmp2 <- geneAnnot[which(geneAnnot$name == iGene)]
                                             tmp3 <- c(GenomicAlignments::seqnames(tmp2),
                                                       GenomicRanges::strand(tmp2)) %>% as.character() %>%
                                                 c(iGene)
                                             tmpDT <- rep(tmp3, txLengths[iGene]) %>%
                                                 matrix(ncol = 3, byrow = T) %>%
                                                 cbind(txtools:::exonBlockGen(iGene, geneAnnot)) %>%
                                                 cbind(seq(1, txLengths[iGene]))
                                             tmpDT <- tmpDT[, c(1, 4, 2, 3, 5)] %>% data.table::data.table() %>%
                                                 magrittr::set_colnames(c("chr", "gencoor",
                                                                          "strand", "gene", "txcoor"))
                                             tmpDT$chr <- as.factor(tmpDT$chr)
                                             tmpDT$gencoor <- as.integer(tmpDT$gencoor)
                                             tmpDT$strand <- as.factor(tmpDT$strand)
                                             tmpDT$gene <- as.factor(tmpDT$gene)
                                             tmpDT$txcoor <- as.integer(tmpDT$txcoor)
                                             return(tmpDT)
                                         }) %>% do.call(what = "rbind")
        txDT[, 1:5] <- genCoorSys
        txDT
    }
    else {
        stop("Genes of txDT are not contained in geneAnnot$name .\n")
    }
}
