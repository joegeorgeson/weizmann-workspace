# This is a quick guide for using [sailor](https://github.com/YeoLab/sailor) on wexac

The intent is to be used as in the paper [Robust single-cell discovery of RNA targets of RNA-binding proteins and ribosomes](https://www.nature.com/articles/s41592-021-01128-0)

## Usage
```
$ cd /home/labs/schwartzlab/Collaboration/programs/sailor
$ ml Singularity
$ singularity exec sailor-1.2.0.img wf_rnaediting2strands.cwl example.json
```

## Initial setup
* `sif` file was downloaded from https://external-collaborator-data.s3.us-west-1.amazonaws.com/public-software/sailor.tar.gz
* git pulled from https://github.com/YeoLab/sailor.git
* In touch with authors via github [here](https://github.com/YeoLab/sailor/issues/15#issuecomment-1448671771)
