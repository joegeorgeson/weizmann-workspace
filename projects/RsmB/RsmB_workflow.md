# This is a workflow to prepare surfmaps for m5C


Step 1: find similar proteins in temp database
```
conda activate diamond
nohup diamond blastp --query /home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/WP_011011787.1.faa --max-target-seqs 0 --threads 12 --db /home/labs/schwartzlab/Collaboration/BrianRoss/GSHCdb/diamondBLAST/GSHT_proteomes.dmnd --outfmt 6 qseqid sseqid pident qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq full_sseq --out /home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/diamond_blastp_out.txt &
```

Step 2: filter in R
```
diamond.hits <-fread("/home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/diamond_blastp_out.txt")
colnames(diamond.hits) <- c("qseqid","sseqid","pident","qlen","slen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qseq","sseq","full_sseq")

diamond.hits.eval.filt <- diamond.hits[diamond.hits$evalue < 1e-10]

seqinr::write.fasta(as.list(diamond.hits.eval.filt$full_sseq),
                    diamond.hits.eval.filt$sseqid,
                    "/home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/colabfold_input.faa")

```

Step 3: Run MMseqs2
```
conda activate my_colabfold
bsub -q molgen-q -n 24,24 -P colabfold_search -R 'span[hosts=1]' -R rusage[mem=12500] colabfold_search --mmseqs /home/labs/schwartzlab/joeg/alphafold/mmseqs/bin/mmseqs /home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/colabfold_input.faa /home/labs/schwartzlab/Collaboration/databases/alphafold/ /home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/colabfold_search/
```
Step 3b: split alignment files into directories of ~100 to pass multiple jobs ...typically jobs of 100 finish on gpu-short with 1 gpu

```
mmseq.dir <- "/home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/colabfold_search/"
colabfold.dir <- "/home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/"

mmseq.files <- list.files(mmseq.dir, pattern=".a3m")

dir.size <- 100
colabfold.cmd <- c()
for(i in 1:ceiling(length(mmseq.files)/dir.size)){
  if(i==1){
    dir.create(paste0(mmseq.dir,"/",i))
    mv.files <- mmseq.files[1:dir.size]
    mv.cmd <- paste0("mv ", mmseq.dir,"/",mv.files, " ", paste0(mmseq.dir,"/",i), "/", mv.files)
    lapply(mv.cmd, function(x) system(x))
  }
  if(i!=1 & i!=ceiling(length(mmseq.files)/dir.size)){
    dir.create(paste0(mmseq.dir,"/",i))
    mv.files <- mmseq.files[((i-1)*dir.size+1):(i*dir.size)]
    mv.cmd <- paste0("mv ", mmseq.dir,"/",mv.files, " ", paste0(mmseq.dir,"/",i), "/", mv.files)
    lapply(mv.cmd, function(x) system(x))
  }
  if(i==ceiling(length(mmseq.files)/dir.size)){
    dir.create(paste0(mmseq.dir,"/",i))
    mv.files <- mmseq.files[((i-1)*dir.size+1):length(mmseq.files)]
    mv.cmd <- paste0("mv ", mmseq.dir,"/",mv.files, " ", paste0(mmseq.dir,"/",i), "/", mv.files)
    lapply(mv.cmd, function(x) system(x))
  }
  
  colabfold.cmd[i] <- paste0("bsub -R rusage[mem=20000] -gpu num=1:j_exclusive=yes -q gpu-long -C 0 ",
                          "-e err_",i," -o out_",i,
                          " colabfold_batch --amber --templates --num-recycle 3 ", mmseq.dir,"/",i, " ",colabfold.dir,"/predictions_",i)
}


fileConn <- file(paste0(colabfold.dir,"/predict_bsub.sh"))
writeLines(c("#!/bin/sh",
             "ml purge",
             "ml restore af",
             "#conda activate my_colabfold",
             colabfold.cmd),
           sep="\n",
           con=fileConn)
close(fileConn)

```

Step 4: Run alphafold by calling `sh predict_bsub.sh` in the previous step 


Step 5: tmAlign relative to Pfu-RsmB

Step 6: run surfmap (locally, singularity via bsub fails for unknown reason)
```
conda activate surfmap
cd /home/labs/schwartzlab/joeg/github/SURFMAP #need relative paths...
singularity run /home/labs/schwartzlab/joeg/sif/surfmap_v1.5.sif -pdb /home/labs/schwartzlab/joeg/alphafold/Pfu_ac4C/predictions_1/1072_relaxed_rank_1_model_2.pdb -proj flamsteed -d /home/labs/schwartzlab/joeg/alphafold/Pfu_ac4C/surfmap/ -tomap electrostatics
```




