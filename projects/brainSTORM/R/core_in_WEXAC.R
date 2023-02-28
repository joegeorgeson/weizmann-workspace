#' Make STAR genome
#'
#' Wrapper for STAR genome generation
#'
#' @param fastaGenome
#' @param bedAnnotation
#' @param outDir
#' @param nCores
#' @param maxReadLength
#'
#' @return
#' @export
#'
#' @examples
mkSTARgenome <- function(fastaGenome, bedAnnotation = NULL, outDir = NULL,
                         nCores = 2, maxReadLength = 101){
    outName <- paste0(strsplit(fastaGenome, split = "/") %>% unlist %>% utils::tail(1),
                      ".STAR")
    if(is.null(outDir)){
        mkTmpDir()
        outDir <- file.path(getwd(), "STORMtmp_dir", outName)
    }else if(!is.null(outDir)){
        outDir <- outDir
    }
    genome <- txtools::tx_load_genome(fastaGenome)
    lGen <- sum(Biostrings::width(genome))
    nRef <- length(Biostrings::width(genome))
    if(nRef > 1000){
        genomChrBinNbits <- floor(min(18,log2(max(lGen/nRef, maxReadLength))))
    }else{
        genomChrBinNbits <- 18
    }
    genomeindexNb <- floor(min(14, log2(lGen)/2 - 1))
    if(is.null(bedAnnotation)){
        com <- paste("module load STAR/2.7.5c &&",
                     "/apps/RH7U2/general/STAR/2.7.5c/bin/Linux_x86_64/STAR",
                     "--runMode genomeGenerate",
                     "--runThreadN", nCores,
                     "--genomeDir", outDir,
                     "--genomeFastaFiles", fastaGenome,
                     "--genomeChrBinNbits", genomChrBinNbits,
                     "--genomeSAindexNbases", genomeindexNb)
        system(com)
    }else{
        tmpGTF <- mkGTF(bedAnnotation = bedAnnotation, outName = tempfile())
        com <- paste("module load STAR/2.7.5c &&",
                     "/apps/RH7U2/general/STAR/2.7.5c/bin/Linux_x86_64/STAR",
                     "--runMode genomeGenerate",
                     "--runThreadN", nCores,
                     "--genomeDir", outDir,
                     "--genomeFastaFiles", fastaGenome,
                     "--sjdbOverhang", maxReadLength - 1,
                     "--sjdbGTFfile", tmpGTF,
                     "--genomeChrBinNbits", genomChrBinNbits,
                     "--genomeSAindexNbases", genomeindexNb)
        system(com)
        invisible(file.remove(tmpGTF))
    }
    outDir
}

#' STAR alignment
#'
#' Wrapper to perform alignment of sequencing reads to a reference genome
#' using STAR (Dobin) and sorting and indexing using Samtools.
#'
#' @param read1Files character. Path to R1 FASTQ files
#' @param STARgenomeDir character. Path to STAR genome to map to
#' @param pairedEnd logical. Indicating if to expect a Read2 FASTQ file
#' File will be automatically looked for but both need to include "_R1", and "_R2"
#' in their respective file names.
#' @param zipped logical. TRUE as default for zipped FASTQ files, will be read with zcat.
#' If set to FALSE FASTQ files will be read with cat, as not zipped.
#' @param nCores numeric. Number of cores used for alignment and sorting processes.
#' @param outFilterMultimapNmax numeric. Threshold for which a read will have multiple
#' mappings. Default is 10. For unique alignment change the value to 1.
#' @param outDir character. Path to output directory
#' @param alignIntronMax numeric. maximum intron size, if 0, max intron size will be determined
#' by default as (2ˆwinBinNbits)*winAnchorDistNbins (See STAR manual for more info)
#' @param alignEndsType character. type of read ends alignment (See STAR manual for more info)
#' @param otherSTARparams character. Additional parameters not covered by the arguments
#' of this function can be added here, separated by spaces.
#' @param dry logical. If set to TRUE the alignment is not performed, only output
#' are the paths to expected output files.
#' @param tmpDir character. Path to directory to be used to generate
#' intermediate files.
#' @param logSumm logical. Set to FALSE so final logs are not appended to
#' *mappingSummary.txt*.
#'
#' @return
#' @export
alignSTAR <- function(read1Files, STARgenomeDir, pairedEnd = TRUE, zipped = TRUE,
                      nCores = 4, alignEndsType = "Local", alignIntronMax = 0,
                      outFilterMultimapNmax = 10, outDir, otherSTARparams = "",
                      dry = FALSE, tmpDir = NULL, logSumm = TRUE){
    if(!all(grepl(patt = "_R1", read1Files))){stop("All read1Files must contain the string '_R1'")}
    rootNames <- lapply(read1Files, function(x){strsplit(x, split = "/") %>% unlist %>%
            tail(1)}) %>% unlist %>% gsub(pattern = "_R1(.)+", replacement = "")
    if(!dry & is.null(tmpDir)){mkTmpDir()}
    if(is.null(tmpDir)){
        tmpDir <- "STORMtmp_dir"
    }
    bamFiles <- file.path(tmpDir, paste0(rootNames, "_Aligned.out.bam"))
    BAM <- c(gsub(bamFiles, pattern = ".bam$", replacement = ".sorted.bam"),
             gsub(bamFiles, pattern = ".bam$", replacement = ".sorted.bam.bai"))
    outBAM <- lapply(gsub(bamFiles, pattern = ".bam$", replacement = ".sorted.bam"),
                     function(x){strsplit(x, split = "/") %>% unlist %>%
                             tail(1)}) %>% unlist %>% file.path(outDir, .)
    if(dry){return(outBAM)}
    if(!dir.exists(outDir)){dir.create(outDir)}
    if(!dir.exists(tmpDir)){dir.create(tmpDir)}
    if(zipped){rFCom <- "zcat"}else if(!zipped){rFCom <- "cat"}
    for(read1F in read1Files){
        if(pairedEnd){
            read2F <- gsub(read1F, pattern = "_R1", replacement = "_R2")
            if(!file.exists(read2F)){stop(read2F, " does not exist.")}
        }else if(!pairedEnd){
            read2F <- ""
        }else{stop("pairedEnd must be logical either TRUE or FALSE")}
        outFPrefix <- strsplit(read1F, split = "/") %>% unlist %>% tail(1) %>%
            gsub(pattern = "R1.fastq.gz", replacement = "") %>%
            gsub(pattern = "R1.fastq", replacement = "")
        if(is.null(tmpDir)){
            outFPrefix <- file.path("STORMtmp_dir", outFPrefix)
        }else{
            outFPrefix <- file.path(tmpDir, outFPrefix)
        }
        com <- paste0("module load STAR/2.7.5c &&",
                      "/apps/RH7U2/general/STAR/2.7.5c/bin/Linux_x86_64/STAR",
                      " --runMode alignReads",
                      " --runThreadN ", nCores,
                      " --genomeDir ", STARgenomeDir,
                      " --readFilesCommand ", rFCom,
                      " --readFilesIn ", read1F, " ", read2F,
                      " --outFileNamePrefix ", outFPrefix,
                      " --outSAMtype ", "BAM Unsorted",
                      " --outFilterMultimapNmax ", outFilterMultimapNmax,
                      " --alignEndsType ", alignEndsType,
                      " --alignIntronMax ", alignIntronMax,
                      " ", otherSTARparams)
        system(com)
    }
    # Alignment Summary Report
    if(logSumm){
        logFiles <- file.path(tmpDir, paste0(rootNames, "_Log.final.out"))
        RES <- lapply(logFiles, function(x){
            utils::read.delim(file = x, header = FALSE, stringsAsFactors = FALSE)
        })
        # Merge in one table
        summary <- lapply(seq_along(RES), function(x){
            RES[[x]][,2]
        }) %>% do.call(what = cbind)
        rownames(summary) <- RES[[1]][,1]
        colnames(summary) <- rootNames
        outReport <- file.path(outDir, "mappingSummary.txt")
        if(file.exists(outReport)){ # Add columns to existing summary report
            tmp <- data.table::fread(outReport, header = TRUE) %>%
                tibble::column_to_rownames("V1")
            utils::write.table(x = cbind(tmp, summary), file = outReport,
                               sep = "\t", quote = F, col.names = NA)
        }else{
            utils::write.table(x = summary, file = outReport, sep = "\t", quote = F,
                               col.names = NA)
        }
    }
    # Sort and index with samtools
    for(file in bamFiles){
        system(paste0("module load samtools/1.9 && ",
                      "/apps/RH7U2/gnu/samtools/1.9/bin/samtools sort -o ",
                      gsub(file, pattern = ".bam$", replacement = ".sorted.bam"),
                      " ", file, " -@", nCores))
        system(paste0("module load samtools/1.9 && ",
                      "/apps/RH7U2/gnu/samtools/1.9/bin/samtools index ",
                      gsub(file, pattern = ".bam$", replacement = ".sorted.bam")))
        system(paste0("rm ", file))
    }
    # Remove all garbage
    if(logSumm){
        garbageSuffix <- c("_SJ.out.tab", "_Log.progress.out", "_Log.out",
                           "_Log.final.out")
    }else{
        garbageSuffix <- c("_SJ.out.tab", "_Log.progress.out", "_Log.out")
    }
    garbage <- file.path(tmpDir, lapply(rootNames, function(x){
        paste0(x, garbageSuffix)
    }) %>% unlist)
    invisible(file.remove(garbage))
    # Move files to output dir
    for(file in BAM){
        system(paste("mv", file, outDir))
    }
    outBAM
}

#' STORM object constructor from META table
#'
#'
#'
#' @param META data.frame Experimental design
#' @param genome
#' @param geneAnnot
#' @param nCores
#' @param verb
#'
#' @return
#' @export
#'
#' @examples
storm_STORM <- function(META, genome = NULL, geneAnnot = NULL, nCores = 1, verb = FALSE){
    if(hasDups(META$id)){
        stop("Not allowed duplicated id variables in META")
    }
    if(verb){cat("Loading RDS files \n")}
    DTL <- parallel::mclapply(mc.cores = nCores, META$RDS, function(x) readRDS(x)) %>%
        magrittr::set_names(META$id)# load RDS files
    # Check if data is in data.table or list format
    if(verb){cat("Checking if data is in data.table or list format \n")}
    if(all(lapply(DTL, function(x) class(x)[1]) %>% unlist() %>% magrittr::equals("data.frame"))){
        DTL <- lapply(DTL, data.table::data.table) %>% magrittr::set_names(META$id)
    }else if(all(lapply(DTL, class) %>% unlist() %>% magrittr::equals("list"))){
        DTL <- lapply(DTL, function(x){
            do.call(x, what = rbind) %>% data.table::data.table()
        }) %>% magrittr::set_names(META$id)
    }
    # Have reference sequence, if not add it
    if(verb){cat("Checking if txDTs have reference sequence \n")}
    reqRefSeq <- lapply(DTL, function(x){"refSeq" %in% names(x)}) %>% unlist %>% magrittr::not()
    if((sum(reqRefSeq) > 0) & (is.null(genome) | is.null(geneAnnot))){
        stop("Data requires reference sequence, genome and geneAnnot arguments must be provided")
    }
    if(sum(reqRefSeq) > 0){
        if(verb){cat("Adding reference sequence... \n")}
        DTL[reqRefSeq] <- parallel::mclapply(mc.cores = nCores, DTL[reqRefSeq], function(DT){
            txtools::tx_split_DT(DT) %>% lapply(function(x){
                txtools::tx_add_refSeqDT(DT = x, genome = genome, geneAnnot = geneAnnot)
            }) %>% txtools::tx_merge_DT()
        })
    }
    # Check uniformity of DTs, if unequal equalize
    if(verb){cat("Checking uniformity of DTs \n")}
    nRowF <- nrow(DTL[[1]])
    identCoors <- parallel::mclapply(mc.cores = nCores, seq_along(DTL)[-1], function(i){
        if(nRowF == nrow(DTL[[i]])){
            all(DTL[[1]][,c("gene", "txcoor")] == DTL[[i]][,c("gene", "txcoor")])
        }else{
            FALSE
        }
    }) %>% unlist() %>% all()
    if(!identCoors){
        if(verb){cat("Unifying txDTs \n")}
        if(is.null(geneAnnot) | is.null(genome)){
            stop("Data needs compatibility adjustment, this requires geneAnnot and
                 genome arguments to be provided")
        }
        DTL <- txtools::tx_unifyTxDTL(DTL, geneAnnot, genome, nCores, type = "union")
    }
    if(!("pos" %in% names(DTL[[1]]))){
        DTL <- lapply(DTL, function(x) txtools::tx_add_pos(x))
    }
    STORM <- list(META = META,
                  DATA = DTL,
                  RES  = NULL)
    if(verb){cat("STORM object generated. \n")}
    return(STORM)
}

#' Confussion matrix plot
#'
#' @param STORM
#' @param title
#' @param pointAlpha
#'
#' @return
#' @export
#'
storm_plot_confMat <- function(STORM, title = "", pointAlpha = 0.5){
    df <- storm_confMatrix4plot(STORM)
    ggplot(df, aes(x = outcome, y = RNAmod)) +
        geom_jitter(alpha = 0.5, width = 0.2, height = 0.2) +
        theme_bw() + ggtitle(title)
}

#' Make bisulphite genome
#'
#' @param fastaGenome character. Path to genome fasta file
#' @param outFile character. Name of output file. Auto adds
#'
#' @return
#' @export
#'
#' @examples
bisGenome <- function(fastaGenome, outFile = "auto"){
    if(outFile == "auto"){
        mkTmpDir()
        fileN <- strsplit(fastaGenome, split = "/") %>% unlist %>% utils::tail(1)
        outName <- paste0(strsplit(fastaGenome, split = "/") %>% unlist %>% utils::tail(1),
                          ".bis")
        outFile <- file.path(getwd(), "STORMtmp_dir", outName)
    }
    gen <- Biostrings::readDNAStringSet(fastaGenome)
    gen_r <- lapply(gen, function(x){
        tmp <- stringr::str_replace_all(x, pattern = "C", replacement = "T")
    })
    seqinr::write.fasta(sequences = gen_r, names = names(gen_r), file.out = outFile, as.string = T)
    outFile
}

#' Make transcriptome
#'
#' Extract transcriptome sequences from genome and generate FASTA file
#'
#' @param fastaGenome
#' @param bedAnnotation
#' @param outFile
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
mkTranscriptome <- function(fastaGenome, bedAnnotation, outFile = "auto", nCores){
    genome <- txtools::tx_load_genome(fastaGenome)
    geneAnnot <- txtools::tx_load_bed(bedAnnotation)
    seqs <- getGeneSeqsfromGenome(genome = genome, geneAnnot = geneAnnot, nCores)
    if(outFile == "auto"){
        mkTmpDir()
        fileN <- strsplit(fastaGenome, split = "/") %>% unlist %>% utils::tail(1)
        outFile <- file.path("STORMtmp_dir", paste0(fileN, ".txOme"))
    }
    names(seqs) <- geneAnnot$name
    Biostrings::writeXStringSet(seqs, filepath = outFile)
    outFile
}

#' Make BED from FASTA Transcriptome
#'
#' Make BED Gene Annotation from a FASTA genome. Considering that the sequences
#' in the FASTA file represent intron-less genes.
#'
#' @param fastaTxOme
#' @param outFile
#'
#' @return
#' @export
#'
#' @examples
mkBedFromFastaTxOme <- function(fastaTxOme, outFile = "auto"){
    fa <- txtools::tx_load_genome(fastaTxOme)
    if(max(Biostrings::width(fa)) > 10000){warning("Maximum sequence length exceeded 10000, make sure sequences represent transcripts.")}
    if(outFile == "auto"){
        mkTmpDir()
        fileN <- strsplit(fastaTxOme, split = "/") %>% unlist %>% utils::tail(1)
        outFile <- file.path("STORMtmp_dir", paste0(fileN, ".bed"))
    }
    tmp <- GenomicRanges::GRanges(seqnames = names(fa),
                                  ranges = IRanges::IRanges(start = 1, width = Biostrings::width(fa)),
                                  strand = "+")
    names(tmp) <- names(fa)
    plyranges::write_bed(tmp, file = outFile)
    outFile
}

#' Library complexity report
#'
#' @param META
#' @param maxExtrapolation
#' @param steps
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
libComplexReport <- function(META, maxExtrapolation = 2.01e6, steps = 1e5, verbose = FALSE){
    if(all(file.exists(META$BAM))){
        for(file in META$BAM){
            com <- paste0("module load libgsl/2.3 && ",
                          "module load preseq/2.0.1 && ",
                          "/apps/RH7U2/gnu/preseq/2.0.1/preseq lc_extrap -P -B ",
                          "-e ", maxExtrapolation, " -s ", steps, " ", file, " -o ",
                          gsub(pattern = ".bam$", replacement = ".lce.txt", x = file),
                          " &")
            system(com)
        }
        Sys.sleep(time = 20) #TODO: Use the package instead
        lastBam <- gsub(pattern = ".bam$", replacement = ".lce.txt", x = META$BAM)
        while(min(difftime(Sys.time(), file.info(lastBam)$mtime, units = "secs"), na.rm = TRUE) < 60){
            Sys.sleep(time = 2)
        }
    }else{
        stop("Files ", paste(META$BAM[!file.exists(META$BAM)], collapse = " "),
             " do not exist")
    }
    if(verbose){cat("DONE: Library complexity reports.")}
}

#' Align using STAR from META
#'
#' @param META
#' @param STARGenome
#' @param STARGenomebis
#' @param pairedEnd
#' @param outDir
#' @param nCores
#' @param tmpDir
#' @param logSumm
#' @param alignIntronMax numeric. maximum intron size, if 0, max intron size will be determined
#' by default as (2ˆwinBinNbits)*winAnchorDistNbins (See STAR manual for more info)
#' @param otherSTARparams character. Additional parameters not covered by the arguments
#' of this function can be added here, separated by spaces.
#' @param force logical. Set to TRUE to realign existing BAM files.
#'
#' @return
#' @export
#'
#' @examples
STARalign_META <- function(META, STARGenome, STARGenomebis, pairedEnd = TRUE,
                           outDir, nCores, tmpDir, logSumm, alignIntronMax = 0,
                           otherSTARparams = "", force = FALSE){
    if(all(file.exists(META$BAM)) & !force){
        warning("Alignment was not performed as all BAM files already ",
                "exist and force was not imposed.")
        return(NULL)
    }else{
        if(!force){
            META <- META[!file.exists(META$BAM), ]
        }
        # STAR Alignment
        if(any(!META$bisAlign)){
            alignSTAR(read1Files = META[META$bisAlign == FALSE, "FASTQ"],
                      pairedEnd = pairedEnd,
                      nCores = nCores,
                      STARgenomeDir = STARGenome,
                      outDir = outDir,
                      tmpDir = tmpDir,
                      logSumm = logSumm,
                      alignIntronMax = alignIntronMax,
                      otherSTARparams = otherSTARparams) %>% invisible()
        }
        if(any(META$bisAlign)){
            alignSTAR(read1Files = META[META$bisAlign == TRUE, "FASTQ"],
                      pairedEnd = pairedEnd,
                      nCores = nCores,
                      STARgenomeDir = STARGenome_bis,
                      outDir = outDir,
                      tmpDir = tmpDir,
                      logSumm = logSumm,
                      alignIntronMax = alignIntronMax,
                      otherSTARparams = otherSTARparams) %>% invisible()
        }
    }
}


# Helpers ######
# Make GTF from BED
mkGTF <- function(bedAnnotation, outName = NULL, source = "user"){
    if(is.null(outName)){
        mkTmpDir()
        tmpDir <- file.path(getwd(), "STORMtmp_dir")
        outName <- file.path(tmpDir, "tmp_geneAnnot.gtf")
    }
    com <- paste0("module load kentUtils/v377 && ",
                  "/apps/RH7U2/general/kentUtils/v377/bin/bedToGenePred ",
                  bedAnnotation, " /dev/stdout | /apps/RH7U2/general/",
                  "kentUtils/v377/bin/genePredToGtf file /dev/stdin ", outName)
    system(com)
    tmp <- utils::read.delim(outName, header = FALSE)
    tmp[,2] <- source
    utils::write.table(tmp, file = outName, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    outName
}
