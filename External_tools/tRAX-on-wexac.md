# This is documentation for using [tRAX](https://github.com/UCSC-LoweLab/tRAX) on wexac

## Run via singulariry
Example 1
```
ml Singularity
singularity run /home/labs/schwartzlab/Collaboration/programs/tRAX/trax_latest.sif trimadapters.py --runname STORM_tRNA --runfile /home/labs/schwartzlab/joeg/Collaborations/Jordan/human_tRNA_28Feb23/tRAX/runfile.txt
singularity run /home/labs/schwartzlab/Collaboration/programs/tRAX/trax_latest.sif processsamples.py --experimentname=STORM --databasename=/home/labs/schwartzlab/Collaboration/programs/tRAX/hg38 --ensemblgtf=/home/labs/schwartzlab/Collaboration/programs/tRAX/genes.gtf --samplefile=/home/labs/schwartzlab/joeg/Collaborations/Jordan/human_tRNA_28Feb23/tRAX/samplefile.txt --exppairs=/home/labs/schwartzlab/joeg/Collaborations/Jordan/human_tRNA_28Feb23/tRAX/samplespairs.txt --nofrag
```

## Initial setup
Download git
```
cd /home/labs/schwartzlab/Collaboration/programs/tRAX
git clone https://github.com/UCSC-LoweLab/tRAX.git
```

Create conda for building databases locally
```
cd /home/labs/schwartzlab/Collaboration/programs/tRAX
conda env create -f trax_env.yaml
```


To pull docker and create `sif`
```
[joeg@epitrans ~]$ ssh docker1
[joeg@docker1 ~]$ docker pull ucsclowelab/trax
[joeg@docker1 ~]$ docker images |grep ucsclowelab/trax
ucsclowelab/trax                latest                            5bd8d0e2a89c   20 months ago       3.37GB
[joeg@docker1 ~]$ docker tag 5bd8d0e2a89c ops:5000/trax
[joeg@docker1 ~]$ docker push ops:5000/trax
[joeg@docker1 ~]$ singularity pull docker://ucsclowelab/trax
[joeg@docker1 ~]$ scp trax_latest.sif joeg@epitrans:/home/labs/schwartzlab/joeg/github/tRAX/
```


build for hg38
```
conda activate trax_env
sh quickdb.bash hg38 .
python maketrnadb.py --databasename=hg38 --genomefile=genome.fa --trnascanfile=hg38-tRNAs-detailed.out --namemap=hg38-tRNAs_name_map.txt
```

