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
