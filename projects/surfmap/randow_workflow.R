library(data.table)
library(ggplot2)

###
### functions
###
{
  order.mats <- function(file.paths, ref.prot){
    this.path <- paste0(surfmap.dir, file.paths)
    names.use <- paste0("WP_",unlist(lapply(unlist(lapply(file.paths, function(x) strsplit(x, "WP")[[1]][3])), function(y) strsplit(y, "_")[[1]][2])))
    these.mats <- lapply(1:length(this.path), function(x) fread(this.path[x]))
    these.mats.2 <-lapply(these.mats, function(x){
      x$svalue <- as.numeric(gsub("Inf", 0, x$svalue))
      return(x)
    })
    ref.mat.index <- grep(int.prot, this.path)
    mat.sum.diff <- unlist(lapply(1:length(this.path), function(x) sum(abs(these.mats.2[[ref.mat.index]]$svalue - these.mats.2[[x]]$svalue), na.rm=T)))
    order.mats <- these.mats[order(mat.sum.diff, decreasing = F)]
    order.mats <- mapply(cbind, order.mats, "prot"=names.use[order(mat.sum.diff, decreasing = F)], SIMPLIFY=F)
    return(order.mats)
  }
  order.diff <- function(file.paths, ref.prot){
    this.path <- paste0(surfmap.dir, file.paths)
    names.use <- paste0("WP_",unlist(lapply(unlist(lapply(file.paths, function(x) strsplit(x, "WP")[[1]][3])), function(y) strsplit(y, "_")[[1]][2])))
    these.mats <- lapply(1:length(this.path), function(x) fread(this.path[x]))
    these.mats.2 <-lapply(these.mats, function(x){
      x$svalue <- as.numeric(gsub("Inf", 0, x$svalue))
      return(x)
    })
    ref.mat.index <- grep(int.prot, this.path)
    mat.sum.diff <- unlist(lapply(1:length(this.path), function(x) sum(abs(these.mats.2[[ref.mat.index]]$svalue - these.mats.2[[x]]$svalue), na.rm=T)))
    order.diff <- mat.sum.diff[order(mat.sum.diff, decreasing = F)]
    names.diff <- names.use[order(mat.sum.diff, decreasing = F)]
    order.dat <- list(order.diff, names.diff)
    names(order.dat) <- c("diff", "prot")
    return(order.dat)
  }
}
###
### variables
###

surfmap.dir <- "/home/labs/schwartzlab/joeg/tmp/sharon/af/surfmap/"
int.prot <- "WP_011011787.1"

###
### number crunching
###
{
surfmap.mats <- list.files(surfmap.dir, pattern="_matrix.txt", recursive = T)

surfmap.mats.elect <- surfmap.mats[grep("electrostatics", surfmap.mats)]
surfmap.mats.ww <- surfmap.mats[grep("wimley_white", surfmap.mats)]
surfmap.mats.kd <- surfmap.mats[grep("kyte_doolittle", surfmap.mats)]
surfmap.mats.stick <- surfmap.mats[grep("stickiness", surfmap.mats)]

elect.dat <- order.mats(surfmap.mats.elect, int.prot)
ww.dat <- order.mats(surfmap.mats.ww, int.prot)
kd.dat <- order.mats(surfmap.mats.kd, int.prot)
stick.dat <- order.mats(surfmap.mats.stick, int.prot)

elect.diff <- order.diff(surfmap.mats.elect, int.prot)
ww.diff <- order.diff(surfmap.mats.ww, int.prot)
kd.diff <- order.diff(surfmap.mats.kd, int.prot)
stick.diff <- order.diff(surfmap.mats.stick, int.prot)

elect.dat.dt <- do.call(rbind, elect.dat)
ww.dat.dt <- do.call(rbind, ww.dat)
kd.dat.dt <- do.call(rbind, kd.dat)
stick.dat.dt <- do.call(rbind, stick.dat)
}
###
### GSHT, TMalign, BLAST intersect
###
{
bac.GSHTdb <- fread("/home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/bacteria_ncbi_temperatures.csv")
arc.GSHTdb <- fread("/home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/archaea_ncbi_temperatures.csv")
euk.GSHTdb <- fread("/home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/eukaryotes_ncbi_temperatures.csv")
colnames(arc.GSHTdb) <- colnames(bac.GSHTdb)
colnames(euk.GSHTdb) <- colnames(bac.GSHTdb)

bac.GSHTdb.filt <- bac.GSHTdb[!is.na(bac.GSHTdb$Temperature)]
arc.GSHTdb.filt <- arc.GSHTdb[!is.na(arc.GSHTdb$Temperature)]
euk.GSHTdb.filt <- euk.GSHTdb[!is.na(euk.GSHTdb$Temperature)]

all.GSHTdb.filt <- rbind(cbind("bacteria", bac.GSHTdb.filt), cbind("archaea", arc.GSHTdb.filt),
                         cbind("eukaryote", euk.GSHTdb.filt))
GSHT.search.names <- unlist(lapply(all.GSHTdb.filt$`Organism Name`, function(x) paste0(strsplit(x, " ")[[1]][1], " ", strsplit(x, " ")[[1]][2])))
all.GSHTdb.filt$names <- GSHT.search.names

## TM blast intersect
TMdat.paths <- list.files("/home/labs/schwartzlab/joeg/tmp/sharon/af/TMalign/", pattern = "TMalign.out", full.names = T)
TMscores <- unlist(lapply(TMdat.paths, function(x) as.numeric(strsplit(readLines(x)[18], " ")[[1]][2])))
TMscores.names <- unlist(lapply(list.files("/home/labs/schwartzlab/joeg/tmp/sharon/af/TMalign/", pattern = "TMalign.out", full.names = F),
                                function(x) strsplit(x, "_r")[[1]][1]))
blast.used <- readAAStringSet("/home/labs/schwartzlab/joeg/tmp/sharon/af/m5C_blasthit.fasta")
blast.all <- fread("/home/labs/schwartzlab/joeg/tmp/sharon/Pfu_m5C_nrDiamond_out_ultraSenstitive_e10.txt")
colnames(blast.all) <- c("qseqid","sseqid","pident","qlen","slen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qseq","sseq","full_sseq","staxids","sscinames","sskingdoms","skingdoms","sphylums","stitle")
TM.match <- cbind(blast.all[match(TMscores.names, blast.all$sseqid)], TMscores)

dim(blast.all)
hist(log(blast.all$pident))

}
###
### order and plot mats (only electrostatics shown)
###
check.range <- 100

{
  range.elect.top <- c(1:check.range)
  range.elect.last <- c(1,(length(elect.dat)-(check.range-2)):length(elect.dat))

  elect.names.top <- elect.diff$prot[range.elect.top]
  elect.diff.top <- elect.diff$diff[range.elect.top]
  elect.dat.top <- elect.dat[range.elect.top]
  elect.dat.dt.top <- do.call(rbind, elect.dat.top)
  elect.dat.dt.top$prot <- factor(elect.dat.dt.top$prot, levels=elect.names.top)
  elect.blast.orgn.top <- TM.match[match(elect.names.top, substr(TM.match$sseqid,1,nchar(TM.match$sseqid)-2))]
  elect.sci.names.top <- unlist(lapply(elect.blast.orgn.top$sscinames, function(x) strsplit(x, ";")[[1]][1]))
  elect.blast.orgn.top$Temp <- all.GSHTdb.filt$Temperature[match(elect.sci.names.top, all.GSHTdb.filt$names)]
  elect.label.top <- paste0(elect.sci.names.top, "_", elect.names.top, "_Temp:", elect.blast.orgn.top$Temp)

  elect.names.last <- elect.diff$prot[range.elect.last]
  elect.diff.last <- elect.diff$diff[range.elect.last]
  elect.dat.last <- elect.dat[range.elect.last]
  elect.dat.dt.last <- do.call(rbind, elect.dat.last)
  elect.dat.dt.last$prot <- factor(elect.dat.dt.last$prot, levels=elect.names.last)
  elect.blast.orgn.last <- TM.match[match(elect.names.last, substr(TM.match$sseqid,1,nchar(TM.match$sseqid)-2))]
  elect.sci.names.last <- unlist(lapply(elect.blast.orgn.last$sscinames, function(x) strsplit(x, ";")[[1]][1]))
  elect.blast.orgn.last$Temp <- all.GSHTdb.filt$Temperature[match(elect.sci.names.last, all.GSHTdb.filt$names)]
  elect.label.last <- paste0(elect.sci.names.last, "_", elect.names.last, "_Temp:", elect.blast.orgn.last$Temp)

  elect.top <- ggplot(data = elect.dat.dt.top, aes(x=absc, y=ord, fill=svalue)) +
    geom_tile() + facet_grid(.~prot) +
    scale_fill_gradientn(colours = c("green","black","red"))
  elect.last <- ggplot(data = elect.dat.dt.last, aes(x=absc, y=ord, fill=svalue)) +
    geom_tile() + facet_grid(.~prot) +
    scale_fill_gradientn(colours = c("green","black","red"))
  cowplot::plot_grid(elect.top, elect.last, nrow=2)
}

dat.out <- cbind(elect.blast.orgn.top$sscinames,elect.blast.orgn.top$sskingdoms,
                 elect.blast.orgn.top$qseqid, elect.blast.orgn.top$sseqid, elect.blast.orgn.top$pident,
                 elect.blast.orgn.top$TMscores, elect.blast.orgn.top$Temp)
colnames(dat.out) <- c("Name", "domain", "query", "hit", "percID", "TMscore", "growthTemp")
dat.out <- data.table(dat.out)
View(dat.out)
###
### distance
###
library(factoextra)

elect.dat.by.row.l <- lapply(elect.dat, function(x) x$svalue)
elect.dat.by.row.l <- lapply(elect.dat.by.row.l, function(x) gsub("Inf", 0, x))
elect.dat.by.row.l <- lapply(elect.dat.by.row.l, function(x) as.numeric(x))
elect.dat.by.row <- data.table(do.call(rbind, elect.dat.by.row.l))
rownames(elect.dat.by.row) <- elect.diff$prot
colnames(elect.dat.by.row) <- paste0("P", 1:dim(elect.dat.by.row)[2])

dim(elect.dat.by.row)

elect.clust <- factoextra::hcut((elect.dat.by.row), k = 3,
                                hc_metric="manhattan", stand = F)
str(elect.clust)
elect.clust$silinfo

factoextra::fviz_dend(elect.clust)
factoextra::fviz_silhouette(elect.clust)
# factoextra::fviz_cluster(elect.clust)

elect.clust
