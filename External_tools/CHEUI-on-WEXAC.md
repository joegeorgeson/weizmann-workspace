# This is a guide for using [CHEUI](https://github.com/comprna/CHEUI) on WEXAC resources

General notes
* For a human transcriptome you'll need ~3TB of disc space **per sample**
* Many steps are sequential but some can be done in parallel
* This will take a solid day of processing time
* The reference 'genome' should be a transcriptome
* My analysis removed rRNA from GENCODE and added the 'Taoka' reference transcripts
* This approach parallelizaes the process by splitting the bam file and recombines files after model1
* Submission to general queue results in pre-emption so submit to schwartz queue ...
* Index fast5 files first as this is rate limiting step

Items can be downloaded and untar-ed as below;
```
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR470/ERR4706156/HEK293T-WT-rep1.tar.gz
ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR470/ERR4706161/HEK293T-WT-rep1.fastq.gz

tar -xzvf HEK293T-WT-rep1.tar.gz
```

I built their conda with cuda inside the conda (the other version failed)
```
conda create --name cheui python=3.7 tensorflow-gpu=2.4.1 pandas=1.3.4 conda-forge::cudatoolkit-dev -y && conda activate cheui
```

The first chunk of analysis is for preparing the fast5 and bam file for processing
```
ml nanopolish/0.13.2
nohup nanopolish index -d /home/labs/schwartzlab/joeg/reviews/CHEUI/aws/ERR4706156/20181120_0819_SHO_20112018_EmptyE2-9/fast5/ /home/labs/schwartzlab/joeg/reviews/CHEUI/aws/ERR4706161/HEK293T-WT-rep1.fastq.gz &

ml minimap2
nohup minimap2 -ax map-ont -k14 --split-prefix -dHEK293T-WT-rep1.tmp -t12 /home/labs/schwartzlab/joeg/genomes/Homo_sapiens_hg38/hg38_rRNA_txOme_exons_good.fa /home/labs/schwartzlab/joeg/reviews/CHEUI/aws/ERR4706161/HEK293T-WT-rep1.fastq.gz -o /home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/HEK293T-WT-rep1.sam &

ml SAMtools/1.15.1-GCC-11.2.0
samtools view -Sb -F 2324 -@ 12 HEK293T-WT-rep1.sam > HEK293T-WT-rep1.2324.bam
samtools sort -@ 12 HEK293T-WT-rep1.2324.bam -o HEK293T-WT-rep1.2324.sorted.bam
samtools index -@ 12 HEK293T-WT-rep1.2324.sorted.bam
rm HEK293T-WT-rep1.sam HEK293T-WT-rep1.2324.bam #cleanup when possible
samtools idxstats HEK293T-WT-rep1.2324.sorted.bam > HEK293T-WT-rep1.2324.sorted.bam.idxstats.txt
```

This next chunk of `R` code was used to generate job submissions;
```
read.path <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/aws/ERR4706161/HEK293T-WT-rep1.fastq.gz"
in.bam <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/HEK293T-WT-rep1.2324.sorted.bam"
genome.path <- "/home/labs/schwartzlab/joeg/genomes/Homo_sapiens_hg38/hg38_rRNA_txOme_exons_good.fa"
working.dir <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/"

read.path <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/aws/keep/HEK293T-WT-rep2.fastq.gz"
in.bam <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep2/HEK293T-WT-rep2.2324.sorted.bam"
genome.path <- "/home/labs/schwartzlab/joeg/genomes/Homo_sapiens_hg38/hg38_rRNA_txOme_exons_good.fa"
working.dir <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep2/"


read.path <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/aws/ERR4706163/HEK293T-WT-rep3.fastq.gz"
in.bam <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep3/HEK293T-WT-rep3.2324.sorted.bam"
genome.path <- "/home/labs/schwartzlab/joeg/genomes/Homo_sapiens_hg38/hg38_rRNA_txOme_exons_good.fa"
working.dir <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep3/"


dir.create(working.dir)
nanopolish.dir <-paste0(working.dir,"/event_align/")
dir.create(nanopolish.dir)
bam.dir <- paste0(working.dir,"/split_bam/")
dir.create(bam.dir)
m6A.dir <- paste0(working.dir,"/m6A_preprocess/")
dir.create(m6A.dir)
cmd.dir <- paste0(working.dir,"/cmds/")
dir.create(cmd.dir)

ff <- fread(paste0(in.bam,".idxstats.txt"))
in.list <- ff$V1[ff$V3>20]


samtools.cmd <- paste0("samtools view -Sb ", in.bam, " ", in.list, " > ", bam.dir, in.list, ".bam")
samtools.idx <- paste0("samtools index ", bam.dir, in.list, ".bam")

bam.files <- paste0(in.list, ".bam")
out.name <- gsub(".bam","_nanopolish_out.txt",bam.files)

eventalign.cmd <- paste0("nanopolish eventalign -t 2 ",
                         "--reads ",read.path, " --bam ", bam.dir, bam.files, " --genome ",
                         genome.path, " --print-read-names --scale-events --samples > ",nanopolish.dir, out.name) # for C++ implementation
# genome.path, " --scale-events --signal-index  --samples --print-read-names") # original implementation


m6A.out.dir <- gsub("_nanopolish_out.txt", "_m6A_preprocess", out.name)

m6A.preprocess.cmd <- paste0("python3 /home/labs/schwartzlab/joeg/github/CHEUI/scripts/CHEUI_preprocess_m6A.py ",
                             "-i ", nanopolish.dir, out.name, 
                             " -o ", m6A.dir, m6A.out.dir,
                             " -m /home/labs/schwartzlab/joeg/github/CHEUI/kmer_models/model_kmer.csv",
                             " -n 2")

m6A.dir.files <- paste0(m6A.dir, m6A.out.dir,"/",in.list, "_nanopolish_out_signals+IDS.p")
m6A.read.out <- gsub("_nanopolish_out.txt","_m6A_read_level.txt",out.name)
local.name <- gsub("_nanopolish_out.txt","",out.name)

m6A.model1.cmd <- paste0("python /home/labs/schwartzlab/joeg/github/CHEUI/scripts/CHEUI_predict_model1.py -i ",
                         m6A.dir.files," -m /home/labs/schwartzlab/joeg/github/CHEUI/CHEUI_trained_models/CHEUI_m6A_model1.h5 ",
                         "-o ",m6A.dir, m6A.read.out, " -l ", local.name)

sort.cmd <- paste0("sort -k1 --parallel=2 ", m6A.dir, m6A.read.out, " > ", m6A.dir, m6A.read.out,".sorted")

m6A.model2.out <- gsub("read_level.txt.sorted", "site_level.txt", paste0(m6A.dir, m6A.read.out,".sorted"))

m6A.model2.cmd <- paste0("python /home/labs/schwartzlab/joeg/github/CHEUI/scripts/CHEUI_predict_model2.py -i ",
                         m6A.dir, m6A.read.out,".sorted -m /home/labs/schwartzlab/joeg/github/CHEUI/CHEUI_trained_models/CHEUI_m6A_model2.h5 ",
                         "-o ", m6A.model2.out)

all.sh.cmds <- paste0(cmd.dir, in.list,".sh")

for(i in 1:length(samtools.cmd)){
  
  fileConn<-file(all.sh.cmds[i])
  writeLines(c("#!/bin/sh",
               samtools.cmd[i],
               samtools.idx[i],
               eventalign.cmd[i],
               m6A.preprocess.cmd[i],
               m6A.model1.cmd[i],
               sort.cmd[i],
               m6A.model2.cmd[i]),
             sep="\n",
             con=fileConn)
  close(fileConn)
}


all.bsub.cmd <- paste0("bsub -q schwartz -P run_CHEUI_rep2 -n 2,2 -R 'span[hosts=1]' -R rusage[mem=500] ",
                       "sh ", all.sh.cmds)

fileConn<-file(paste0(working.dir,"/all_cmds_bsub.sh"))
writeLines(c("#!/bin/sh",
             "ml purge",
             "ml nanopolish/0.13.2",
             "ml SAMtools/1.15.1-GCC-11.2.0",
             "#conda activate cheui2",
             all.bsub.cmd),
           sep="\n",
           con=fileConn)
close(fileConn)


print("run in terminal;")
print("conda activate cheui2")
print(paste0("nohup ", working.dir, "/all_cmds_bsub.sh &"))
```

lastly the data can be consolidated as below;
```
HEK293.WT.rep1.m6A <- get.m6A.site("/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/m6A_preprocess/")
HEK293.WT.rep2.m6A <- get.m6A.site("/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep2/m6A_preprocess/")
HEK293.WT.rep3.m6A <- get.m6A.site("/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep3/m6A_preprocess/")


get.m6A.site <- function(in.dir){
  site.files <- list.files(in.dir, pattern="site", full.names = T)
  site.dat <- lapply(site.files, function(x) fread(x))
  site.dat.all <- rbindlist(site.dat)
  return(site.dat.all)
}
```
