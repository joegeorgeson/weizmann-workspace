# This is a guide for using [CHEUI](https://github.com/comprna/CHEUI) on WEXAC resources

General notes
* For a human transcriptome you'll need ~3TB of disc space **per sample**
* Many steps are sequential but some can be done in parallel
* This will take a solid day of processing time
* The reference genome should be a transcriptome
* My analysis removed rRNA from GENCODE and added the 'Taoka' reference transcripts
* This approach parallelizaes the process by splitting the bam file and recombines files after model1

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

The first chunk of analysis is for preparing bam files and indexing the fast5 files
```
ml minimap2
bsub -q new-short -P minimap2 -n 24,24 -R 'span[hosts=1]' -R rusage[mem=10000] minimap2 -ax map-ont -k14 --split-prefix -dHEK293T-WT-rep1tmp -t24 /home/labs/schwartzlab/joeg/genomes/Homo_sapiens_hg38/hg38_rRNA_txOme_exons_good.fa /home/labs/schwartzlab/joeg/reviews/CHEUI/aws/ERR4706161/HEK293T-WT-rep1.fastq.gz -o /home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/HEK293T-WT-rep1.sam

ml SAMtools/1.15.1-GCC-11.2.0
samtools view -Sb -F 2324 -@ 12 HEK293T-WT-rep1.sam > HEK293T-WT-rep1.2324.bam
samtools sort -@ 12 HEK293T-WT-rep1.2324.bam -o HEK293T-WT-rep1.2324.sorted.bam
samtools index -@ 12 HEK293T-WT-rep1.2324.sorted.bam

rm HEK293T-WT-rep1.sam HEK293T-WT-rep1.2324.bam #cleanup when possible

ml picard/2.25.1-Java-11
mkdir split_bam
java -jar $EBROOTPICARD/picard.jar SplitSamByNumberOfReads -I HEK293T-WT-rep1.2324.sorted.bam -O split_bam -N_FILES 200

ml SAMtools/1.15.1-GCC-11.2.0
cd split_bam
bamf=*bam
for f in $bamf; do samtools index -@12 $f; done

ml nanopolish/0.13.2
nanopolish index -d /home/labs/schwartzlab/joeg/reviews/CHEUI/aws/ERR4706156/20181120_0819_SHO_20112018_EmptyE2-9/fast5/ /home/labs/schwartzlab/joeg/reviews/CHEUI/aws/ERR4706161/HEK293T-WT-rep1.fastq.gz #note there was no sequencing summary in the aws files, will take few hours
```

This next chunk of `R` code was used to generate job submissions;
```
nanopolish.dir <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/event_align/"
read.path <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/aws/ERR4706161/HEK293T-WT-rep1.fastq.gz"
genome.path <- "/home/labs/schwartzlab/joeg/genomes/Homo_sapiens_hg38/hg38_rRNA_txOme_exons_good.fa"

bam.dir <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/split_bam/"
bam.files <- list.files(bam.dir, pattern="bam$", full.names=F)
out.name <- gsub(".bam","_nanopolish_out.txt",bam.files)

eventalign.cmd <- paste0("bsub -q new-short -P nanopolish_eventalign -n 12,12 -R 'span[hosts=1]' -R rusage[mem=1000] -N -oo ",
                         nanopolish.dir, out.name, " nanopolish eventalign -t 12 ",
                         "--reads ",read.path, " --bam ", bam.dir, bam.files, " --genome ",
                         genome.path, " --print-read-names --scale-events --samples") # for C++ implementation
# genome.path, " --scale-events --signal-index  --samples --print-read-names") # original implementation

fileConn<-file("/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/nanopolish_eventalign_bsub.sh")
writeLines(c("#!/bin/sh",
             "ml purge",
             "ml nanopolish/0.13.2",
             eventalign.cmd),
           sep="\n",
           con=fileConn)
close(fileConn)
```

Now we need to use the 'preprocessing' scripts from CHEUI for m6A and m5C, again from R. Note that the m6A worked better with the python script and the m5C worked better with the C++ script
```
m6A.dir <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/m6A_preprocess/"
m6A.out.dir <- gsub("_nanopolish_out.txt", "_m6A_preprocess", out.name)

m5C.dir <- "/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1/m5C_preprocess/"
m5C.out.dir <- gsub("_nanopolish_out.txt", "_m5C_preprocess", out.name)

m6A.preprocess.cmd <- paste0("bsub -q new-short -P m6A_preprocess -n 8,8 -R 'span[hosts=1]' -R rusage[mem=1000] ",
                             "'python3 /home/labs/schwartzlab/joeg/github/CHEUI/scripts/CHEUI_preprocess_m6A.py ",
                             "-i ", nanopolish.dir, out.name, 
                             " -o ", m6A.dir, m6A.out.dir,
                             " -m /home/labs/schwartzlab/joeg/github/CHEUI/kmer_models/model_kmer.csv",
                             " -n 8'")

m5C.preprocess.cmd <- paste0("bsub -q new-short -P m5C_preprocess -n 8,8 -R 'span[hosts=1]' -R rusage[mem=1000] ",
                             "'/home/labs/schwartzlab/joeg/github/CHEUI/scripts/preprocessing_CPP/CHEUI ",
                             "-i ", nanopolish.dir, out.name, 
                             " -o ", m5C.dir, m5C.out.dir,
                             " -m /home/labs/schwartzlab/joeg/github/CHEUI/kmer_models/model_kmer.csv",
                             " -n 8 --m5C'")

fileConn<-file("/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1//preprocess_m6A_bsub.sh")
writeLines(c("#!/bin/sh",
             "ml purge",
             "#conda activate cheui2",
             m5C.preprocess.cmd),
           sep="\n",
           con=fileConn)
close(fileConn)

fileConn<-file("/home/labs/schwartzlab/joeg/reviews/CHEUI/HEK293T-WT-rep1//preprocess_m5C_bsub.sh")
writeLines(c("#!/bin/sh",
             "ml purge",
             "ml SAMtools/1.15.1-GCC-11.2.0", #this is the version of C I used ...it worked
             "#conda activate cheui2",
             "#cd /home/labs/schwartzlab/joeg/github/CHEUI/scripts/preprocessing_CPP",
             m5C.preprocess.cmd),
           sep="\n",
           con=fileConn)
close(fileConn)
```

to be continued...


