Downloaded the references from `http://ftp.ebi.ac.uk/pub/databases/alphafold/`

```
wget http://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.csv
wget http://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta
```


Download specific pdb files
```
wget https://alphafold.ebi.ac.uk/files/{alphafold_ID}-model_{database_version}.pdb

# for example;
wget https://alphafold.ebi.ac.uk/files/AF-A0A009GC07-F1-model_v4.pdb
```
