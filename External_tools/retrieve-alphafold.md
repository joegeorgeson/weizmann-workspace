# This is a guide for retreieving alphafold files from the ftp site



## Download specific files
```
# pdb
wget https://alphafold.ebi.ac.uk/files/{alphafold_ID}-model_{database_version}.pdb
# cif
wget https://alphafold.ebi.ac.uk/files/{alphafold_ID}-model_{database_version}.cif
# bcif
wget https://alphafold.ebi.ac.uk/files/{alphafold_ID}-model_{database_version}.bcif
# bcif
wget https://alphafold.ebi.ac.uk/files/{alphafold_ID}-model_{database_version}.bcif
# pae json
wget https://alphafold.ebi.ac.uk/files/{alphafold_ID}-predicted_aligned_error_{database_version}.bcif
# pae png
wget https://alphafold.ebi.ac.uk/files/{alphafold_ID}-predicted_aligned_error_{database_version}.png
```

## for example;
```
wget https://alphafold.ebi.ac.uk/files/AF-A0A009GC07-F1-model_v4.pdb
wget https://alphafold.ebi.ac.uk/files/AF-A0A009GC07-F1-model_v4.cif
wget https://alphafold.ebi.ac.uk/files/AF-A0A009GC07-F1-model_v4.bcif
wget https://alphafold.ebi.ac.uk/files/AF-A0A009GC07-F1-predicted_aligned_error_v4.json
wget https://alphafold.ebi.ac.uk/files/AF-A0A009GC07-F1-predicted_aligned_error_v4.png
```


## Seach accession IDs
This file can be loaded into R and accessions matched with Uniprot
```
/home/labs/schwartzlab/joeg/alphafold/diamond/accession_ids.rds
```

## You can do a blast search (with dimaond) against the below
```
/home/labs/schwartzlab/joeg/alphafold/diamond/alphafold.dmnd
```


## Downloaded the references from http://ftp.ebi.ac.uk/pub/databases/alphafold/

```
wget http://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.csv
wget http://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta
```
