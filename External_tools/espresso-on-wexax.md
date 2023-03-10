# Below are guidelines for [ESPRESSO](https://github.com/Xinglab/espresso) on wexac

## About
ESPRESSO (Error Statistics PRomoted Evaluator of Splice Site Options) is a novel method for processing alignment of long read RNA-seq data, which can effectively improve splice junction accuracy and isoform quantification.

## Installation
* This is to run locally (without snakemake)
* Straightforward to install via conda requirements
```
cd /home/labs/schwartzlab/Collaboration/programs
git clone https://github.com/Xinglab/espresso.git
cd /home/labs/schwartzlab/Collaboration/programs/espresso/snakemake
conda create --name espresso -c conda-forge bioconda --file conda_requirements.txt
```

## Usage
* The below shows an example of the test samples via the command line (not snakemake)
* All step are to be done inside the conda `conda activate espresso`

### Step 1: prepare bed and align fastq
```
conda activate espresso
cd /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44
paftools.js gff2bed cd44.gtf > cd44.junctions.bed
minimap2 -ax splice -ub -k14 -w 4 --junc-bed cd44.junctions.bed --secondary=no cd44.fasta PC3E_1_cd44.fastq > PC3E_1_cd44.sam
minimap2 -ax splice -ub -k14 -w 4 --junc-bed cd44.junctions.bed --secondary=no cd44.fasta PC3E_2_cd44.fastq > PC3E_2_cd44.sam
minimap2 -ax splice -ub -k14 -w 4 --junc-bed cd44.junctions.bed --secondary=no cd44.fasta PC3E_3_cd44.fastq > PC3E_3_cd44.sam
minimap2 -ax splice -ub -k14 -w 4 --junc-bed cd44.junctions.bed --secondary=no cd44.fasta GS689_1_cd44.fastq > GS689_1_cd44.sam
minimap2 -ax splice -ub -k14 -w 4 --junc-bed cd44.junctions.bed --secondary=no cd44.fasta GS689_2_cd44.fastq > GS689_2_cd44.sam
minimap2 -ax splice -ub -k14 -w 4 --junc-bed cd44.junctions.bed --secondary=no cd44.fasta GS689_3_cd44.fastq > GS689_3_cd44.sam
samf=*sam
for f in $samf; do samtools sort -o $f.sorted.sam $f; done
```

### Step 2: create samples.tsv like the below
```
/home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/PC3E_1_cd44.sam.sorted.sam	PC3E_1_cd44
/home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/PC3E_2_cd44.sam.sorted.sam	PC3E_2_cd44
/home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/PC3E_3_cd44.sam.sorted.sam	PC3E_3_cd44
/home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/GS689_1_cd44.sam.sorted.sam	GS689_1_cd44
/home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/GS689_2_cd44.sam.sorted.sam	GS689_2_cd44
/home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/GS689_3_cd44.sam.sorted.sam	GS689_3_cd44
```

### Step 3: process via espresso
```
# Step 1
perl /home/labs/schwartzlab/Collaboration/programs/espresso/src/ESPRESSO_S.pl -L /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/samples.tsv -F /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.fasta -A /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.gtf -O /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir

# Step 2
perl /home/labs/schwartzlab/Collaboration/programs/espresso/src/ESPRESSO_C.pl -I /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/ -F /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.fasta -X 0
perl /home/labs/schwartzlab/Collaboration/programs/espresso/src/ESPRESSO_C.pl -I /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/ -F /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.fasta -X 1
perl /home/labs/schwartzlab/Collaboration/programs/espresso/src/ESPRESSO_C.pl -I /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/ -F /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.fasta -X 2
perl /home/labs/schwartzlab/Collaboration/programs/espresso/src/ESPRESSO_C.pl -I /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/ -F /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.fasta -X 3
perl /home/labs/schwartzlab/Collaboration/programs/espresso/src/ESPRESSO_C.pl -I /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/ -F /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.fasta -X 4
perl /home/labs/schwartzlab/Collaboration/programs/espresso/src/ESPRESSO_C.pl -I /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/ -F /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.fasta -X 5

# Step 3
perl /home/labs/schwartzlab/Collaboration/programs/espresso/src/ESPRESSO_Q.pl -L /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/samples.tsv.updated -A /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.gtf
```

### Step 4: Visualize
```
python3 /home/labs/schwartzlab/Collaboration/programs/espresso/visualization/visualize.py --genome-fasta /home/labs/schwartzlab/Collaboration/programs/espresso/test_data/test_data_espresso_cd44/cd44.fasta --updated-gtf /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/samples_N2_R0_updated.gtf --abundance-esp /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/samples_N2_R0_abundance.esp --target-gene cd44 --minimum-count 1 --descriptive-name cd44 --output-dir /home/labs/schwartzlab/Collaboration/programs/espresso/work_dir/visualization
```


