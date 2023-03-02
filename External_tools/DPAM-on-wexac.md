This is a guide for using [DPAM](https://github.com/CongLabCode/DPAM) on wexac


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

The conda doesn't correctly source eveything, so this is the workaround for now.
```
export PATH=/home/labs/schwartzlab/Collaboration/programs/DPAM/mkdssp:/home/labs/schwartzlab/Collaboration/programs/DPAM/DaliLite.v5/bin/dali.pl:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/bin:/home/labs/schwartzlab/Collaboration/programs/DPAM/hh-suite/build/scripts:/home/labs/schwartzlab/joeg/.conda/envs/DPAM/bin:$PATH
```
