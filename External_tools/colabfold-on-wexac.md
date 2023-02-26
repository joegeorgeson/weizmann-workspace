# [colabfold](https://github.com/sokrypton/ColabFold) was adopted from their github page

## Step 0: Build their database
Use their `setup_databases.sh` to download and setup everything making sure to use the correct version of `mmseqs` and on a machine mimicking the memory/cpu used in `colabfold_search` so the index sizes are compatible.

Downloaded and setup to (1.5TB): `/home/labs/schwartzlab/Collaboration/databases/alphafold`

## Step 1: MSA of all proteins with MMseqs2 (make sure to use their binary)
```
conda activate my_colabfold
bsub -q molgen-q -n 24,24 -P colabfold_search -R 'span[hosts=1]' -R rusage[mem=12500] colabfold_search --mmseqs /home/labs/schwartzlab/joeg/alphafold/mmseqs/bin/mmseqs /home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/colabfold_input.faa /home/labs/schwartzlab/Collaboration/databases/alphafold/ /home/labs/schwartzlab/joeg/alphafold/Pfu_m5C/colabfold_search/
```

## Step 2: Organize job submissions and fold on gpu nodes
```
conda activate my_colabfold



