# This is some documentation for analyzing STAMP data for Sharon

## Retrieve genome & annotation for ensembl
* Used as a convenience for functions, can port over to UCSC/NCBI if needed
```
wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
```

## Create transcriptome (easier to have reference outside of R)
```
ml BEDTools/2.30.0-GCC-11.2.0
bedtools getfasta -split -fi Homo_sapiens.GRCh37.dna.primary_assembly.fa -bed Homo_sapiens.GRCh37.75.gtf -s -name > Homo_sapiens.GRCh37.dna.primary_assembly.tx.fa
```

