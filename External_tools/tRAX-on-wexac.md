
```
git clone https://github.com/UCSC-LoweLab/tRAX.git
```

Create environment
```


build for hg38
```
conda activate trax_env
sh quickdb.bash hg38 .
python maketrnadb.py --databasename=hg38 --genomefile=genome.fa --trnascanfile=hg38-tRNAs-detailed.out --namemap=hg38-tRNAs_name_map.txt
```

more to come
