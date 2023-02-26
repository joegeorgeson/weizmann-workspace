library(data.table)
library(parallel)
### GHSTdb
bac.GSHTdb <- fread("/home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/bacteria_ncbi_temperatures.csv")
arc.GSHTdb <- fread("/home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/archaea_ncbi_temperatures.csv")
euk.GSHTdb <- fread("/home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/eukaryotes_ncbi_temperatures.csv")
colnames(arc.GSHTdb) <- colnames(bac.GSHTdb)
colnames(euk.GSHTdb) <- colnames(bac.GSHTdb)

bac.GSHTdb.filt <- bac.GSHTdb[!is.na(bac.GSHTdb$Temperature)]
arc.GSHTdb.filt <- arc.GSHTdb[!is.na(arc.GSHTdb$Temperature)]
euk.GSHTdb.filt <- euk.GSHTdb[!is.na(euk.GSHTdb$Temperature)]


all.GSHTdb.filt <- rbind(cbind("bac", bac.GSHTdb.filt), cbind("arc", arc.GSHTdb.filt), cbind("euk", euk.GSHTdb.filt))
colnames(all.GSHTdb.filt)[1] <- "domain"

bac.ass <- unlist(lapply(bac.GSHTdb.filt$`GenBank FTP`, function(x) strsplit(x, "GCA_")[[1]][2]))
arc.ass <- unlist(lapply(arc.GSHTdb.filt$`GenBank FTP`, function(x) strsplit(x, "GCA_")[[1]][2]))
euk.ass <- unlist(lapply(euk.GSHTdb.filt$`GenBank FTP`, function(x) strsplit(x, "GCA_")[[1]][2]))

bac.cmd <- paste0("wget ",bac.GSHTdb.filt$`GenBank FTP`, "/GCA_", bac.ass,"_genomic.fna.gz")
arc.cmd <- paste0("wget ",arc.GSHTdb.filt$`GenBank FTP`, "/GCA_", arc.ass,"_genomic.fna.gz")
euk.cmd <- paste0("wget ",euk.GSHTdb.filt$`GenBank FTP`, "/GCA_", euk.ass,"_genomic.fna.gz")

setwd("/home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/bac_genomic/")
mclapply(bac.cmd, function(x) system(x), mc.cores=12)
setwd("/home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/arc_genomic/")
mclapply(arc.cmd, function(x) system(x), mc.cores=12)
setwd("/home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/euk_genomic/")
mclapply(euk.cmd, function(x) system(x), mc.cores=12)
