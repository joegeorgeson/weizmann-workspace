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

## Just download from alphafold!
```
aa <- readRDS("/home/labs/schwartzlab/joeg/alphafold/accession_ids.rds")


diamond.hits <-fread("/home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/diamond_blastp_out.txt")
colnames(diamond.hits) <- c("qseqid","sseqid","pident","qlen","slen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qseq","sseq","full_sseq")

diamond.hits.uniprot.id <- unlist(lapply(diamond.hits$sseqid, function(x) strsplit(x, "_")[[1]][3]))

write.csv(diamond.hits.uniprot.id, row.names=F, "/home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/diamond_blastp_out_names.txt")

uniprot.hits <- fread("/home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.02.27-20.49.45.88.tsv")

af.hits <- aa[aa$Uniprot.ID %in% uniprot.hits$Entry]

wget.cmd <- paste0("wget https://alphafold.ebi.ac.uk/files/",af.hits$AF.ID,"-model_v", af.hits$version, ".pdb")

fileConn<-file("/home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/wget_pdb_cmds.sh")
writeLines(c("#!/bin/sh",
             wget.cmd,
           sep="\n"),
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



### IN PROGRESS
```
pfu.m5C <- "1072_relaxed_rank_1"
pdb.in.dir <- "/home/labs/schwartzlab/joeg/alphafold/Pfu_ac4C/predictions_1/"
mat.in.dir <- "/home/labs/schwartzlab/joeg/alphafold/Pfu_ac4C/TMalign/"


pdb.files <- list.files(pdb.in.dir, pattern="pdb$", full.names = F, recursive = F)
pdb.files <- pdb.files[grep("relaxed_rank_1",pdb.files)]
pfu.m5C.file <- pdb.files[grep(pfu.m5C, pdb.files)]

tmalign.py <- readLines("/home/labs/schwartzlab/joeg/scripts/tmalign_source.py")

load.mobile <- paste0('cmd.load("',pdb.in.dir,pdb.files,'")')
load.target <- paste0('cmd.load("',pdb.in.dir,pfu.m5C.file,'")')
align.cmd <- paste0('tmalign("',
                    substr(pdb.files,1,nchar(pdb.files)-4),'", "WP_011011787.1_relaxed_rank_1_model_3")')
save.cmd <- paste0('cmd.save("',
                   mat.in.dir, substr(pdb.files,1,nchar(pdb.files)-4),'_align.pdb",',
                   '"',substr(pdb.files,1,nchar(pdb.files)-4),'",',0,')')

### this writes the pymol (python) scripts
for(i in 1:length(load.mobile)){
  fileConn<-file(paste0("//home/labs/schwartzlab/joeg/alphafold/Pfu_ac4C/cmds/align_cmd_",i,".py")) #update the path
  writeLines(c(tmalign.py,
               load.target,
               load.mobile[i],
               align.cmd[i],
               save.cmd[i]),
             sep="\n",
             con=fileConn)
  close(fileConn)
}



tm.py.bsub <- paste0("bsub -q new-short -P pymol_tmalign -n 1,1 -R 'span[hosts=1]' -R rusage[mem=1000] ",
                     "pymol -c /home/labs/schwartzlab/joeg/tmp/sharon/af/all_input_PDB/pymol/cmds/align_cmd_",1:1386,".py")

### this is the bsub command to execute the pymol scripts
fileConn<-file("/home/labs/schwartzlab/joeg/tmp/sharon/af/all_input_PDB/pymol/cmds/bsub_align_cmd.sh")
writeLines(c("#!/bin/sh",
             "ml pymol/1.8.6", # this
             tm.py.bsub),
           sep="\n",
           con=fileConn)
close(fileConn)


pdb.files.aligned <- list.files("/home/labs/schwartzlab/joeg/tmp/sharon/af/all_input_PDB/pymol/", pattern="pdb")
bsub.singularity.surfmap.all.cmd <- paste0("bsub -q new-short -P surfmap_singularity -n 1,1 -R 'span[hosts=1]' -R rusage[mem=1000] ",
                                           "singularity run /home/labs/schwartzlab/joeg/sif/surfmap_v1.5.sif -pdb /home/labs/schwartzlab/joeg/tmp/sharon/af/all_input_PDB/pymol/",
                                           pdb.files.aligned,
                                           " -tomap all")

fileConn<-file("/home/labs/schwartzlab/joeg/tmp/sharon/af/all_input_PDB/pymol/cmds/bsub_surfmap_singularity.sh")
writeLines(c("#!/bin/sh",
             bsub.singularity.surfmap.all.cmd),
           sep="\n",
           con=fileConn)
close(fileConn)

singularity.surfmap.all.cmd <- paste0("singularity run /home/labs/schwartzlab/joeg/sif/surfmap_v1.5.sif -pdb /home/labs/schwartzlab/joeg/tmp/sharon/af/all_input_PDB/pymol/",
                                      pdb.files.aligned,
                                      " -tomap electrostatics")

fileConn<-file("/home/labs/schwartzlab/joeg/tmp/sharon/af/all_input_PDB/pymol/cmds/surfmap_singularity_elect.sh")
writeLines(c("#!/bin/sh",
             singularity.surfmap.all.cmd),
           sep="\n",
           con=fileConn)
close(fileConn)
```



