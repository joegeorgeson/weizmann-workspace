#### updated 11-May-2023

# General NGS workflow for Weizmann sequencing

Move raw data from INCPM to local directory via `wget`
```
wget -l 7 -nH --cut-dirs=1 -r --reject='index.html*' --no-parent --no-check-certificate https://stefan.weizmann.ac.il/raw_nov_incpm/PATH_TO_RUN/
```

An example of our routine `bcl2fastq` script is below
* For NovaSeq runs it's typically necessary to apply the mask to remove the 1st base from the second index
  * Using the command as shown, I2 will need to have the first base removed in the SampleSheet.csv
* All indexes have at least 2 base differences in  I1 and I2 so allowing 1 mismatch is fine
* This should take <20 minutes for a full 800M SP100 kit, typically setup as 61|8|8|61 (but can vary depending on requirements of library)


```
ml bcl2fastq/2.20.0.422
bsub -q new-short -P bcl2fastq -n 24,24 -R 'span[hosts=1]' -R rusage[mem=1000] bcl2fastq --runfolder-dir /home/labs/schwartzlab/joeg/data/lib708/230101_A00929_0841_AHNY2CDRX2/ -p 24 --output-dir /home/labs/schwartzlab/joeg/data/lib708/230101_A00929_0841_AHNY2CDRX2/parental_fastq_mask/ --sample-sheet /home/labs/schwartzlab/joeg/data/lib708/SampleSheet_lib708.csv --use-bases-mask Y*,I*,n1i7,Y* --no-lane-splitting --barcode-mismatches 1,1
```


Below is a simple command to evaluate the top 10 indexes
```
gzip -cd gzip_FASTQ_files.fastq.gz | grep '^@' | cut -d: -f10 | sort | uniq -c | sort -k1rn | head -n10
```
