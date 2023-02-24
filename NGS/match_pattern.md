vin.R1 <- readFastq("/home/labs/schwartzlab/vinithi/Projects/Yeast/Engineering_modifications_on_rRNA/2Ome/Seq_Ult_Feb2023/BC10_G2_S77_R1_001.fastq.gz")
vin.R2 <- readFastq("/home/labs/schwartzlab/vinithi/Projects/Yeast/Engineering_modifications_on_rRNA/2Ome/Seq_Ult_Feb2023/BC10_G2_S77_R2_001.fastq.gz")

R1.pattern <- "TGCTTTTGCTGGTAAAAGAAAAACTCTTGTTTCAAT"

R1.match <- vmatchPattern(pattern = R1.pattern, subject = sread(vin.R1), max.mismatch = 1)
index.match <- unlist(lapply(R1.match@ends, function(x) is.null(x)==F))

vin.R1.match <- vin.R1[index.match]
vin.R2.match <- vin.R2[index.match]

vin.R1.NOTmatch <- vin.R1[!index.match]
vin.R2.NOTmatch <- vin.R2[!index.match]

ShortRead::writeFastq(vin.R1.match, "/home/labs/schwartzlab/joeg/tmp/vinithra/14Feb/R1_matches.fastq.gz")
ShortRead::writeFastq(vin.R2.match, "/home/labs/schwartzlab/joeg/tmp/vinithra/14Feb/R2_matches.fastq.gz")

ShortRead::writeFastq(vin.R1.NOTmatch, "/home/labs/schwartzlab/joeg/tmp/vinithra/14Feb/R1_NOTmatches.fastq.gz")
ShortRead::writeFastq(vin.R2.NOTmatch, "/home/labs/schwartzlab/joeg/tmp/vinithra/14Feb/R2_NOTmatches.fastq.gz")
