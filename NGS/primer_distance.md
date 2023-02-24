The below was created in `R` to find the least similar sequences to the existing set of m5C barcodes.

* One approach maximizes the complexity of the sequence based on acss::entropy and acss::entropy2
* The second approach maximizes the hamming distance based on DescTools::StrDist

Recall that the first set of m5C barcodes in the Schwartz lab were mistakingly ordered with the 'NNN' after the barcode instead of before. This still allows for a UMI but may result in some adapters preferentially ligating based on sequence (still not clear if this actually happens.

Some thought should be placed on this because modified RNA oligos cost ~$4k each

```
library(acss)
library(DescTools)
library(data.table)
library(Biostrings)

num.add <- 4
seq.len <- 7
GC.cutoff <- 0.3
# base.space <- c("A","C","T","G")
base.space <- c("A","T","G")
existing.seqs <- c("AATTAGT","AGAATGG","AGTATTG","ATGAAAG","AATGTGA","AGTGATA","AATGGGT","AGGTAGA")
sort.priority <- "complexity"
sort.priority <- "hamming.distance"

new.primers <- primer.dist(num.add, seq.len, base.space, existing.seqs, GC.cutoff,sort.priority="hamming.distance")
new.primers

new.primers <- primer.dist(num.add, seq.len, base.space, existing.seqs, GC.cutoff,sort.priority="complexity")
new.primers


primer.dist <- function(num.add, seq.len, base.space, existing.seqs,GC.cutoff,sort.priority="complexity"){
  base.space.l <- lapply(1:seq.len, function(x) x <- base.space)
  all.seqs <- do.call(paste0,expand.grid(base.space.l))
  
  all.seqs.trip <- regexpr("([[:alpha:]])\\1{1,}",ignore.case = T, all.seqs,perl=T) # finds duplicate repeats
  all.seqs <- all.seqs[attr(all.seqs.trip,"capture.start") == -1] # removes duplicate repeats
  
  all.seq.SS <- DNAStringSet(all.seqs)
  
  GC.content <- letterFrequency(all.seq.SS,letters="GC")/width(all.seq.SS)
  all.seq.SS <- all.seq.SS[as.numeric(GC.content[,1]) >= GC.cutoff]
  
  for(i in 1:(num.add+1)){
    if(i == 1 ){
      seqs.take <- c()
      seqs.take <- existing.seqs
    }else{
      
      seq.distances <- lapply(as.character(all.seq.SS), function(x) lapply(seqs.take, function(y) DescTools::StrDist(x,y, method="hamming")))
      seq.distances.l <- unlist(do.call(rbind, seq.distances))
      seq.distances.dt <- data.table(matrix(seq.distances.l, ncol=length(seqs.take)))
      colnames(seq.distances.dt) <- seqs.take
      seq.distances.dt$totalDist <- rowSums(seq.distances.dt)
      seq.distances.dt$new.seqs <- as.character(all.seq.SS)
      seq.distances.dt.sort <- seq.distances.dt[order(seq.distances.dt$totalDist, decreasing=T)]
      seq.distances.dt.sort$suitable <- apply(seq.distances.dt.sort,1,function(x) all(x>3))
      seq.distances.dt.sort$seq.complexity.1 <- acss::entropy(seq.distances.dt$new.seqs)
      seq.distances.dt.sort$seq.complexity.2 <- acss::entropy2(seq.distances.dt$new.seqs)
      seq.distances.dt.sort$seq.complexity.total <- as.numeric(seq.distances.dt.sort$seq.complexity.1) + as.numeric(seq.distances.dt.sort$seq.complexity.2)
      
      seq.distances.dt.sort <- seq.distances.dt.sort[seq.distances.dt.sort$suitable == T]
      
      if(sort.priority=="complexity"){
        seq.distances.dt.sort <- seq.distances.dt.sort[with(seq.distances.dt.sort, order(-seq.complexity.total, -totalDist))]
      }
      if(sort.priority=="hamming.distance"){
        seq.distances.dt.sort <- seq.distances.dt.sort[with(seq.distances.dt.sort, order(-totalDist, -seq.complexity.total))]
      }
      
      seqs.take <- c(seqs.take, seq.distances.dt.sort$new.seqs[1])
    }
    
  }
  
  return.dat <- list(seqs.take[!seqs.take %in% existing.seqs], existing.seqs)
  names(return.dat) <- c("new.seqs", "existing.seqs")
  return(return.dat)
  
}


###
### from terminal
###
{
  seq.reps <- paste0(rep("{A,T,C,G}",seq.len),collapse="")
  seq.cmd <- paste0("echo ",seq.reps, " > all.seqs.txt")
  seq.cmd <- paste0("echo ",seq.reps)
  
  # system(seq.cmd, intern=T)
  # run seq.cmd from terminal ...system(seq.cmd, intern=T) doesn't work and not troubleshooting
  # all.seqs <- fread("/home/labs/schwartzlab/joeg/tmp/m5C_oligos/all.seqs.txt")
  # all.seq.SS <- DNAStringSet(names(all.seqs))
  
}
```
