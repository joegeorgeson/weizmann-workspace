
```
git clone https://github.com/UCSC-LoweLab/tRAX.git
```

Create environment
```
```

To pull docker;
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

To run singularity;
```
ml Singularity
singularity run /home/labs/schwartzlab/Collaboration/programs/tRAX/trax_latest.sif trimadapters.py --runname STORM_tRNA --runfile /home/labs/schwartzlab/joeg/Collaborations/Jordan/human_tRNA_28Feb23/tRAX/runfile.txt
```
