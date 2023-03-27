# Documentation for analyzing STAMP data for Sharon

## Analyze data
* The below shows an example output from the Rscript
```
library(ggseqlogo)
STAMP.data <- readRDS("/home/labs/schwartzlab/joeg/tmp/sharon/STAMP/joeSTAMP_ex.rds")
STAMP.data <- STAMP.data[as.numeric(STAMP.data$mismatch) / as.numeric(STAMP.data$depth) > 0.1]

STAMP.data <- STAMP.data[STAMP.data$num.hits > 1]

STAMP.seqs.reps <- rep(STAMP.data$seq, as.numeric(STAMP.data$depth))

c.dna <- unlist(lapply(STAMP.seqs.reps, function(x) strsplit(x, ",")[[1]]))
c.dna <- c.dna[nchar(c.dna) == names(which(table(nchar(c.dna)) == max(table(nchar(c.dna)))))]
p1 <- ggseqlogo(c.dna, method = 'bits' )
p2 <- ggseqlogo(c.dna, method = 'prob' )
gridExtra::grid.arrange(p1, p2)
```

## Setup
### Retrieve genome & annotation for ensembl
* Used as a convenience for ensembl functions, can port over to UCSC/NCBI or other GRCh37/hg19 annotation as needed
```
wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
```

### Create transcriptome (easier to have reference outside of R)
* Created bed from gtf with ucsc tools before creating transcriptome below
```
ml BEDTools/2.30.0-GCC-11.2.0
bedtools getfasta -split -fi Homo_sapiens.GRCh37.dna.primary_assembly.fa -bed Homo_sapiens.GRCh37.75.gtf.genePred.bed -s -name > Homo_sapiens.GRCh37.dna.primary_assembly.tx.fa
```

### Run the below in R
* Relies on the pre-built above transcriptome
* Makes use of the `ensembldb::genomeToTranscript` function to retrieve mRNA coordinates and associated transcripts for each genomic coordinate

```
library(ensembldb)
library(EnsDb.Hsapiens.v75)
library(data.table)
library(ggseqlogo)

tx.fa.path <- "/home/labs/schwartzlab/joeg/genomes/Homo_sapiens_ensembl_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.tx.fa"
bed.path.p <- "/home/labs/schwartzlab/joeg/tmp/sharon/STAMP/RBFOX2_STAMP_high_trimmed_refbase_alignUnmapped_final_filterAligned.sortedByCoord.out.rev.sorted.rmdup.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed"
bed.path.n <- "/home/labs/schwartzlab/joeg/tmp/sharon/STAMP/RBFOX2_STAMP_high_trimmed_refbase_alignUnmapped_final_filterAligned.sortedByCoord.out.fwd.sorted.rmdup.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed"

bed.data.p <- fread(bed.path.p, header = F)
bed.data.p$V6 <- "+"
bed.data.n <- fread(bed.path.n, header = F)
bed.data.n$V6 <- "-"
all.bed.data <- rbind(bed.data.p, bed.data.n)
all.bed.data <- all.bed.data[all.bed.data$V1 %in% c(paste0("chr",1:22),"chrX", "chrY", "chrM")]
colnames(all.bed.data) <- c("chr", "start", "end", "info", "score", "strand")

all.bed.data$depth <- lapply(all.bed.data$info, function(x) strsplit(x, "[|]")[[1]][1])
all.bed.data$conversion <- lapply(all.bed.data$info, function(x) strsplit(x, "[|]")[[1]][2])
all.bed.data$mismatch <- as.numeric(all.bed.data$depth)*as.numeric(lapply(all.bed.data$info, function(x) strsplit(x, "[|]")[[1]][3]))

##
## filter the data based on some basic criteria
##
all.bed.data <- all.bed.data[all.bed.data$depth > 20 & all.bed.data$mismatch > 4]

tx.fa.dna <- readDNAStringSet(tx.fa.path)
names(tx.fa.dna) <- unlist(lapply(names(tx.fa.dna), function(x) strsplit(x, "::")[[1]][1]))

hg19.hit.s <- split(all.bed.data, all.bed.data$chr)

window.size <- 25

for(i in 1:length(hg19.hit.s)){
  
  hg19.hit.s[[i]]$seq <- ""
  hg19.hit.s[[i]]$tx.name <- ""
  hg19.hit.s[[i]]$num.hits <- ""
  hg19.hit.s[[i]]$mRNA.pos <- ""
  
  this.bed.dat <- hg19.hit.s[[i]] 
  
  this.bed.dat <- this.bed.dat[order(this.bed.dat$chr)]
  
  hg19.dbx <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == names(hg19.hit.s)[i])
  seqlevelsStyle(hg19.dbx) <- "UCSC"
  
  this.bed.dat.gr <- GRanges(seqnames = as.character(this.bed.dat$chr),
                             ranges=as.character(this.bed.dat$start+1),
                             strand = as.character(this.bed.dat$strand))
  
  tx.gr <- ensembldb::genomeToTranscript(this.bed.dat.gr, hg19.dbx)
  
  
  for(x in 1:length(tx.gr)){
    
    this.gr <- tx.gr[[x]]
    index.match <- which(paste0(hg19.hit.s[[i]]$end,":",hg19.hit.s[[i]]$strand) == unique(paste0(this.gr@elementMetadata@listData$seq_start, ":", this.gr@elementMetadata@listData$seq_strand)))
    
    if(any(this.gr@start==-1)){next}
    if(length(grep("ENST", names(this.gr)))<1){next}
    all.window.seq <- c()
    for(z in 1:length(this.gr)){
      this.window.seq <- c()
      if(length(grep("ENST", names(this.gr)[z]))<1){next}
      this.tx.dna <- tx.fa.dna[grep(names(this.gr)[z],names(tx.fa.dna))]
      win.end <- this.gr@start[z]+window.size
      if(win.end > width(this.tx.dna)){win.end <- width(this.tx.dna)}
      win.start <- this.gr@start[z]-window.size
      if(this.gr@start[z]-window.size < 1){win.start <- 1}
      this.window.seq <- subseq(this.tx.dna, start=win.start, end=win.end)
      all.window.seq <- c(all.window.seq, this.window.seq)
      
    }
    all.window.seq <- do.call(c, all.window.seq)
    
    hg19.hit.s[[i]]$num.hits[index.match] <- length(all.window.seq)
    hg19.hit.s[[i]]$mRNA.pos[index.match] <- paste0(this.gr@start, collapse=",")
    hg19.hit.s[[i]]$seq[index.match] <- paste0(as.character(all.window.seq),collapse = ",")
    hg19.hit.s[[i]]$tx.name[index.match] <- paste0(names(all.window.seq),collapse = ",")
  }
  saveRDS(hg19.hit.s[[i]], file=paste0("/home/labs/schwartzlab/joeg/tmp/sharon/STAMP/", names(hg19.hit.s)[i], "_joeSTAMP.rds"))
}
```
