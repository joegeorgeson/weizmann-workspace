# This is guidance for using [mim-tRNAseq](https://github.com/nedialkova-lab/mim-tRNAseq) on wexac resources

At time of creation, this requires running both in a conda and loading some additional modules. Maybe I'll package into `sif` in the future.


## Installation
* Their git has documentation which should be checked, but at the time or writing this page the conda installer doesn't install everything required for use on wexac.
* Create a conda from the `yaml` in the directory

```
cd /home/labs/schwartzlab/Collaboration/programs/mim-tRNAseq
conda env create -f mimseq_wexac.yml
cp usearch ~/.conda/envs/mimseq/bin/ #needed to anyway source from usearch module ...didn't investigate further
```


## Usage

* Make sure to update the `--out-dir` path
* In Oct-23 it was noted that the current Modomics data this package retrieves is not compatible and users should add the flag `--local-modomics`
* Also in Oct-23 updated the conda env to `mimseq_yaml` because of issues with their conda/github install
* Below is a working example at the time of writing this document (output in folder `hg38_HEK239vsK562`)

```
cd /home/labs/schwartzlab/Collaboration/programs/mim-tRNAseq
ml purge
conda activate mimseq
module load HTSlib/1.12-GCC-10.2.0 samtools/1.9 ncbi-blast+ infernal bedtools/2.26.0 usearch/10.0 gmap/2019-02-26

# Example command
mimseq --species Hsap --cluster-id 0.97 --threads 15 --min-cov 0.0005 --max-mismatches 0.075 --control-condition HEK293T -n hg38_test --out-dir UPDATE_ME --max-multi 4 --remap --remap-mismatches 0.05 sampleData_HEKvsK562.txt
```




