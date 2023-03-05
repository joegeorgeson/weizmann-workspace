# This is a guide for using [DPAM](https://github.com/CongLabCode/DPAM) on wexac


## Environment setup
It will be easiest to setup your own conda. The last step creates a symlink I found necessary.
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
conda install -c salilab dssp
conda install -c anaconda libboost=1.67.0
ln -s /home/labs/schwartzlab/Collaboration/programs/DPAM/DaliLite.v5/bin/dali.pl ~/.conda/envs/DPAM/bin/
ln -s /home/labs/schwartzlab/Collaboration/programs/DPAM/mkdssp ~/.conda/envs/DPAM/bin/
```
And lastly, the conda doesn't correctly source eveything, so this is the workaround for now.
```
export PATH=/home/labs/schwartzlab/Collaboration/programs/DPAM/mkdssp:/home/labs/schwartzlab/Collaboration/programs/DPAM/DaliLite.v5/bin/dali.pl:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/bin:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/scripts:$PATH
```

Check installation
```
$ conda activate DPAM
$ export PATH=/home/labs/schwartzlab/Collaboration/programs/DPAM/mkdssp:/home/labs/schwartzlab/Collaboration/programs/DPAM/DaliLite.v5/bin/dali.pl:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/bin:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/scripts:$PATH
$ cd /home/labs/schwartzlab/Collaboration/programs/DPAM
$ python check_dependencies.py
HH-suite, Foldseek and dali.pl are found
```

## Usage
Here is the general command
```
python DPAM.py <input_cif/pdb> <input_pae> <accession> <output_dir> <threads> <datadir>
```

And something that worked for me (make sure to use json and not png)
```
conda activate DPAM
export PATH=/home/labs/schwartzlab/Collaboration/programs/DPAM/mkdssp:/home/labs/schwartzlab/Collaboration/programs/DPAM/DaliLite.v5/bin/dali.pl:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/bin:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/scripts:$PATH
wget https://alphafold.ebi.ac.uk/files/AF-Q9ZFH0-F1-model_v4.pdb
wget https://alphafold.ebi.ac.uk/files/AF-Q9ZFH0-F1-predicted_aligned_error_v4.json
python /home/labs/schwartzlab/Collaboration/programs/DPAM/DPAM.py AF-Q9ZFH0-F1-model_v4.pdb AF-Q9ZFH0-F1-predicted_aligned_error_v4.json Q9ZFH0 ./ 12 /home/labs/schwartzlab/Collaboration/databases/DPAM/
```


## Database setup (only done once)
They provide a script `download_all_data.sh` to download everything, but I just did it in parallel as below and untar with `tar -xzvf *tgz`
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


