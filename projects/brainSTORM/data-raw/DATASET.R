## code to create brainSTORM's in-package objects

# Other objects ################################################################

storm_baseCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos",
                    "refSeq", "set", "nuc")

storm_baseCoorCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos",
                        "refSeq")

RNAmods_vec <- c("Y", "Am", "Cm", "Gm", "Tm", "m5C", "ac4C", "m1A", "m7G",
                 "m3U", "m62A", "m1acp3Y")

usethis::use_data(storm_baseCols, overwrite = TRUE)
usethis::use_data(storm_baseCoorCols, overwrite = TRUE)
usethis::use_data(RNAmods_vec, overwrite = TRUE)

# RNAmod METRICS ######################################################################
# Metrics to be assigned to specific RNAmods
Y_metrics <- c("SRD1bpDS_CMC.TGIRT_Mock.TGIRT",
               "SRD1bpDS_CMC.SSIII_Mock.SSIII",
               "SRlog2FCh1bpDS_CMC.TGIRT_Mock.TGIRT",
               "SRlog2FCh1bpDS_CMC.SSIII_Mock.SSIII",
               "DRD_m5C.TGIRT_Mock.TGIRT",
               "DRD_RBSseqHeatMg.SSIII_Mock.SSIII",
               "DRD_RBSseqHeatMg.TGIRT_Mock.TGIRT")

Nm_metrics <- c("NmStopScore_MocklowdNTPs.SSIII_Mock.SSIII",
                "ScoreA3p_Mock.TGIRT",
                "ScoreA3p_Mock.SSIII",
                "ScoreA3p_Mock.RTHIV",
                "ScoreC3p_Mock.RTHIV",
                "ScoreC3p_Mock.SSIII",
                "ScoreC3p_Mock.TGIRT")

ac4C_metrics <- c("CtoT.MRD_Ac4C.TGIRT_Mock.TGIRT",
                  "CtoT.MRD_Ac4C.TGIRT_DeacetylatedAc4C.TGIRT",
                  "CtoT.MRD_Ac4C.SSIII_Mock.SSIII",
                  "CtoT.MRD_Ac4C.SSIII_DeacetylatedAc4C.SSIII",
                  "CytPer_m5C.TGIRT_Mock.TGIRT") # Metric suggested by analysis

m1A_metrics <- c("SRD1bpDS_Mock.SSIII_Dimroth.SSIII",
                 "SRlog2FCh1bpDS_Mock.SSIII_Dimroth.SSIII",
                 "MRD_Mock.TGIRT_Dimroth.TGIRT",
                 "MRD_Mock.SSIII_RBSseqHeatMg.SSIII",
                 "MRD_Mock.TGIRT_RBSseqHeatMg.TGIRT",
                 "MRD_Mock.TGIRT_AlkBmix.TGIRT",
                 "ScoreC3p_Mock.RTHIV", # Metric suggested by analysis
                 "ScoreC3p_Mock.TGIRT") # Metric suggested by analysis

m7G_metrics <- c("MRD_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                 "MRD_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                 "SRD1bpDS_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                 "SRD1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                 "SRlog2FCh1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                 "SRlog2FCh1bpDS_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                 "MRD_NaBH4HydBiotin.RTHIV_Mock.RTHIV",
                 "MRD_NaBH4HydBiotin.TGIRT_Mock.TGIRT")

m5C_metrics <- c("CytPer_m5C.TGIRT_Mock.TGIRT",
                 "CytPer_RBSseqHeatMg.SSIII_Mock.SSIII",
                 "CytPer_RBSseqHeatMg.TGIRT_Mock.TGIRT")

m3U_metrics <- c("MRD_Mock.TGIRT_AlkBmix.TGIRT",
                 "ScoreC3p_Mock.TGIRT") # Metric suggested by analysis

m1acp3Y_metrics <- c("SRD1bpDS_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                     "SRD1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                     "SRD1bpDS_Mock.SSIII_Dimroth.SSIII",
                     "SRD1bpDS_CMC.SSIII_Mock.SSIII")

m62A_metrics <- c("MRD_Mock.SSIII_RBSseqHeatMg.SSIII", # Metric suggested by analysis
                  "MRD_Mock.TGIRT_RBSseqHeatMg.TGIRT", # Metric suggested by analysis
                  "SRlog2FCh1bpDS_Mock.SSIII_Dimroth.SSIII", # Metric suggested by analysis
                  "DRD_RBSseqHeatMg.TGIRT_Mock.TGIRT", # Metric suggested by analysis
                  "MRDAtoG_Mock.TGIRT_m5C.TGIRT", # Metric suggested by analysis
                  "MRDAtoG_Mock.SSIII_m5C.SSIII") # Metric suggested by analysis

Y_scores <- Y_metrics
Nm_scores <- Nm_metrics
ac4C_scores <- ac4C_metrics
m1A_scores <- m1A_metrics
m7G_scores <- m7G_metrics
m5C_scores <- m5C_metrics
m3U_scores <- m3U_metrics
usethis::use_data(Y_scores, overwrite = TRUE)
usethis::use_data(Nm_scores, overwrite = TRUE)
usethis::use_data(ac4C_scores, overwrite = TRUE)
usethis::use_data(m1A_scores, overwrite = TRUE)
usethis::use_data(m7G_scores, overwrite = TRUE)
usethis::use_data(m5C_scores, overwrite = TRUE)
usethis::use_data(m3U_scores, overwrite = TRUE)

usethis::use_data(Y_metrics, overwrite = TRUE)
usethis::use_data(Nm_metrics, overwrite = TRUE)
usethis::use_data(ac4C_metrics, overwrite = TRUE)
usethis::use_data(m1A_metrics, overwrite = TRUE)
usethis::use_data(m7G_metrics, overwrite = TRUE)
usethis::use_data(m5C_metrics, overwrite = TRUE)
usethis::use_data(m3U_metrics, overwrite = TRUE)
usethis::use_data(m1acp3Y_metrics, overwrite = TRUE)
usethis::use_data(m62A_metrics, overwrite = TRUE)

# Metrics table
RNAmod_metrics <- lapply(names(metricsList), function(RNAmod){
    cbind(RNAmod, metric = metricsList[[RNAmod]])
}) %>% do.call(what = "rbind") %>% data.table::data.table()
RNAmod_metrics$baseNuc <- unlist(RNAMod_baseNuc[RNAmod_metrics$RNAmod])
RNAmod_metrics$metricFunction <- unlist(base::strsplit(RNAmod_metrics$metric, "_") %>% lapply(function(x) x[1]))
RNAmod_metrics$RTase <- stringr::str_extract(string = RNAmod_metrics$metric, pattern = "(TGIRT|SSIII|RTHIV)")
RNAmod_metrics$groupA <- unlist(base::strsplit(RNAmod_metrics$metric, "_") %>% lapply(function(x) x[2])) %>%
    stringr::str_extract(pattern = paste0("(", paste(vList$libTreat, collapse = "|"), ")"))
RNAmod_metrics$groupB <- unlist(base::strsplit(RNAmod_metrics$metric, "_") %>% lapply(function(x) x[3])) %>%
    stringr::str_extract(pattern = paste0("(", paste(vList$libTreat, collapse = "|"), ")"))
RNAmod_metrics <- RNAmod_metrics[order(RNAmod_metrics$baseNuc), ]
colnames(RNAmod_metrics)
RNAmod_metrics <- RNAmod_metrics[,c(1, 3, 2, 6, 7, 5, 4)]
usethis::use_data(RNAmod_metrics, overwrite = TRUE)

allMetrics <- unique(unlist(metricsList))
usethis::use_data(allMetrics, overwrite = TRUE)
# META default variables list ##################################################

vList <- list(libNum = c("499", "524", "553"),
              organism = c("PyroAbyss", "TherAcid", "Yeast", "Human", "Hs_Sc", "ThermoKoda"),
              RTase = c("SSIII", "SSIV", "TGIRT", "RTHIV"),
              libTreat = c("Ac4C", "CMC", "DeacetylatedAc4C", "MocklowdNTPs",
                           "Mock", "Dimroth", "m5C", "AlkBmix", "RBSseqHeatMg",
                           "NaBH4HydBiotin"),
              bioTreat = c("65deg", "75deg", "85deg", "95deg", "80deg", "95deg", "100deg", "AcidpH1", "AcidpH2", "AcidpH3"),
              replicate = c("rep1", "rep2", "rep3"))
usethis::use_data(vList, overwrite = TRUE)

# Make sample STORM object. ####################################################

# Genome and Gene Annotation
fastaGenome <- "/home/labs/schwartzlab/miguelg/BIGDATA/genome_references/ribosomal/Sc_ribosomal_seqs.fa"
bedAnnotation <- "/home/labs/schwartzlab/miguelg/BIGDATA/gene_annotations/ribosomal/Sc_ribosomal_seqs.bed"
r1Files <- grep(list.files("/home/labs/schwartzlab/miguelg/WORKSPACE_wexac/storm_seq/smpDATA",
                           full.names = TRUE), pattern = "R1", value = TRUE)
OUTDIR <- "/home/labs/schwartzlab/miguelg/WORKSPACE_wexac/storm_seq/lib524/s_cerevisiae"
EXP_NAME <- "S.cerevisiae_SK1-lib524"
NCORES <- 10

# Experimental design table
META <- storm_META(fileNames = r1Files,
                   varsList = vList,
                   outDir = OUTDIR)

# Make transcriptome, make new gene annotation for transcriptome
fastaTxOme <- mkTranscriptome(fastaGenome, bedAnnotation, nCores = NCORES, outFile = "data/fastaTxOme")
bedTxOme <- mkBedFromFastaTxOme(fastaTxOme, outFile = "data/bedTxOme")

# Creating bisulphite transcriptome
bisTxPath <- bisGenome(fastaTxOme, outFile = "data/bisTxPath")

# Create STAR genomes
STARGenome <- mkSTARgenome(fastaTxOme, bedTxOme, outDir = "data/STARGenome")
STARGenome_bis <- mkSTARgenome(bisTxPath, bedTxOme, outDir = "data/STARGenome_bis")

#Loading transcriptome and annotation
GENOME <- txtools::tx_load_genome(fastaTxOme)
TXOME <- txtools::tx_load_bed(bedTxOme)

if(!all(file.exists(META$BAM))){
    # STAR Alignment
    alignSTAR(read1Files = META[libTreat != "m5C" & libTreat != "RBSseqHeatMg", FASTQ],
              nCores = NCORES,
              zipped = TRUE,
              STARgenomeDir = STARGenome,
              alignEndsType = "Local",
              outSAMtype = "BAM Unsorted",
              outDir = OUTDIR)
    alignSTAR(read1Files = META[libTreat == "m5C" | libTreat == "RBSseqHeatMg", FASTQ],
              nCores = NCORES,
              zipped = TRUE,
              STARgenomeDir = STARGenome_bis,
              alignEndsType = "Local",
              outSAMtype = "BAM Unsorted",
              outDir = OUTDIR)
}

if(!all(file.exists(META$RDS))){
    rdsFiles <- parallel::mclapply(mc.cores = NCORES, seq_along(META$FASTQ), function(i){
        bam2TxDT(BAMfile = META$BAM[i],
                 geneAnnot = TXOME,
                 genome = GENOME,
                 dtType = "covNuc",
                 outDir = OUTDIR,
                 nCores = 1,
                 remL = 1000,
                 minR = 0)
    })
}
yeast_STORM <- storm_STORM(META, GENOME, TXOME, nCores = 1)
# Save object to .rda file
usethis::use_data(yeast_STORM, overwrite = TRUE)
# Remove tmp files and dirs
tmpDirs <- c(STARGenome, STARGenome_bis)
tmpFiles <- c(fastaTxOme, bisTxPath, bedTxOme)
rmdFiles <- file.remove(tmpFiles)
unlink(c(STARGenome, STARGenome_bis), recursive = TRUE)

# Taoka Sc sites ###############################################################
rRNAmods_Sc_Taoka <- readRDS("/home/labs/schwartzlab/miguelg/BIGDATA/RNAmod_Annot/Taoka/rib_mods_Sc.rds")
usethis::use_data(rRNAmods_Sc_Taoka, overwrite = TRUE)

# Internal objects #############################################################

# Metrics list
metricsList <- list(Y = Y_metrics, Am = Nm_metrics, Cm = Nm_metrics, Gm = Nm_metrics,
                    Tm = Nm_metrics, ac4C = ac4C_metrics, m1A = m1A_metrics,
                    m7G = m7G_metrics, m5C = m5C_metrics, m3U = m3U_metrics,
                    m62A = m62A_metrics, m1acp3Y = m1acp3Y_metrics)

# RNAMod functions list
RNAModFunList <- list(Y = is.pseudoU, Am = is.Am, Cm = is.Cm, Gm = is.Gm,
                         Tm = is.Tm, m5C = is.m5C, ac4C = is.ac4C, m1A = is.m1A,
                         m7G = is.m7G, m3U = is.m3U, m62A = is.m62A,
                         m1acp3Y = is.m1acp3Y, m6A = is.m6A)

RNAMod_baseNuc <- list(Y = "T", Am = "A", Cm = "C", Gm = "G", Tm = "T",
                       ac4C = "C", m1A = "A", m7G = "G", m5C = "C",
                       m3U = "T", Nm = c("A", "C", "G", "T"), m62A = "A",
                       m1acp3Y = "T", m6A = "A")

# Need to update all internals with one call of use_data(internal = TRUE)!!!
use_data(RNAModFunList, metricsList, RNAMod_baseNuc, internal = TRUE, overwrite = TRUE)
