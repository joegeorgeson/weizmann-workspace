# This is a guide for using [DPAM](https://github.com/CongLabCode/DPAM) on wexac


## Environment setup
It will be easiest to setup your own conda
```
conda create --name DPAM python=3.8
conda activate DPAM
conda install -c qianlabcode dpam
conda install -c bioconda hhsuite
conda install -c ostrokach-forge pdbx
conda install -c bioconda foldseek
conda install -c conda-forge mmcif_pdbx
conda install -c bioconda blast
pip install dssp-wsl
```
And lastly, the conda doesn't correctly source eveything, so this is the workaround for now.
```
export PATH=/home/labs/schwartzlab/Collaboration/programs/DPAM/mkdssp:/home/labs/schwartzlab/Collaboration/programs/DPAM/DaliLite.v5/bin/dali.pl:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/bin:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/scripts:/home/labs/schwartzlab/joeg/.conda/envs/DPAM/bin:$PATH
```

Also, I needed to create a symlink for `dali.pl` in the conda
```
ln -s /home/labs/schwartzlab/Collaboration/programs/DPAM/DaliLite.v5/bin/dali.pl /home/labs/schwartzlab/joeg/.conda/envs/DPAM/bin/
```


## Database setup
They provide a script `download_all_data.sh` to download everything, but I just did it in parallel as below
```
cd /home/labs/schwartzlab/Collaboration/databases/DPAM
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/ECOD70.tgz &
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/pdb70.tgz &
nohup wget --no-check-certificate https://wwwuser.gwdg.de/~compbiol/uniclust/2022_02/UniRef30_2022_02_hhsuite.tar.gz &
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/ECOD_foldseek_DB.tgz &
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/ecod_weights.tgz &
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/ecod_domain_info.tgz &
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/ECOD_benchmark.tgz &
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/ECOD_norms &
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/ecod.latest.domains &
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/ECOD_length &
nohup wget --no-check-certificate https://conglab.swmed.edu/DPAM/ECOD_pdbmap &
```

