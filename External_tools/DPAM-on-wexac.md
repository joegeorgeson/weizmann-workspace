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

And lastly a yml for the conda
```
name: DPAM
channels:
  - qianlabcode
  - anaconda
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - aria2=1.36.0=h1e4e653_3
  - blast=2.13.0=hf3cf87c_0
  - bzip2=1.0.8=h7f98852_4
  - c-ares=1.18.1=h7f98852_0
  - ca-certificates=2023.01.10=h06a4308_0
  - curl=7.88.1=hdc1c0ab_0
  - dpam=2.0=1
  - entrez-direct=16.2=he881be0_1
  - foldseek=5.53465f0=pl5321hf1761c0_0
  - gawk=5.1.0=h7f98852_0
  - gettext=0.21.1=h27087fc_0
  - hhsuite=3.3.0=py38pl5321h8ded8fe_5
  - icu=58.2=he6710b0_3
  - keyutils=1.6.1=h166bdaf_0
  - krb5=1.20.1=h81ceb04_0
  - ld_impl_linux-64=2.40=h41732ed_0
  - libboost=1.67.0=h46d08c1_4
  - libcurl=7.88.1=hdc1c0ab_0
  - libedit=3.1.20191231=he28a2e2_2
  - libev=4.33=h516909a_1
  - libffi=3.4.2=h7f98852_5
  - libgcc-ng=12.2.0=h65d4601_19
  - libgomp=12.2.0=h65d4601_19
  - libiconv=1.17=h166bdaf_0
  - libidn2=2.3.4=h166bdaf_0
  - libnghttp2=1.52.0=h61bc06f_0
  - libnsl=2.0.0=h7f98852_0
  - libsqlite=3.40.0=h753d276_0
  - libssh2=1.10.0=hf14f497_3
  - libstdcxx-ng=12.2.0=h46fd767_19
  - libunistring=0.9.10=h7f98852_0
  - libuuid=2.32.1=h7f98852_1000
  - libxml2=2.9.14=h74e7548_0
  - libzlib=1.2.13=h166bdaf_4
  - lz4-c=1.9.4=h6a678d5_0
  - mmcif_pdbx=2.0.0=pyh44b312d_0
  - ncurses=6.3=h27087fc_1
  - openssl=3.0.8=h0b41bf4_0
  - pcre=8.45=h9c3ff4c_0
  - perl=5.32.1=2_h7f98852_perl5
  - perl-archive-tar=2.40=pl5321hdfd78af_0
  - perl-carp=1.38=pl5321hdfd78af_4
  - perl-common-sense=3.75=pl5321hdfd78af_0
  - perl-compress-raw-bzip2=2.201=pl5321h87f3376_1
  - perl-compress-raw-zlib=2.105=pl5321h87f3376_0
  - perl-encode=3.19=pl5321hec16e2b_1
  - perl-exporter=5.72=pl5321hdfd78af_2
  - perl-exporter-tiny=1.002002=pl5321hdfd78af_0
  - perl-extutils-makemaker=7.66=pl5321hd8ed1ab_0
  - perl-io-compress=2.201=pl5321h87f3376_0
  - perl-io-zlib=1.14=pl5321hdfd78af_0
  - perl-json=4.10=pl5321hdfd78af_0
  - perl-json-xs=2.34=pl5321h9f5acd7_5
  - perl-list-moreutils=0.430=pl5321hdfd78af_0
  - perl-list-moreutils-xs=0.430=pl5321hec16e2b_1
  - perl-parent=0.236=pl5321hdfd78af_2
  - perl-pathtools=3.75=pl5321hec16e2b_3
  - perl-scalar-list-utils=1.62=pl5321hec16e2b_1
  - perl-types-serialiser=1.01=pl5321hdfd78af_0
  - pip=23.0.1=pyhd8ed1ab_0
  - python=3.8.16=he550d4f_1_cpython
  - python_abi=3.8=3_cp38
  - readline=8.1.2=h0f457ee_0
  - setuptools=67.4.0=pyhd8ed1ab_0
  - tk=8.6.12=h27826a3_0
  - wget=1.20.3=ha35d2d1_1
  - wheel=0.38.4=pyhd8ed1ab_0
  - xz=5.2.6=h166bdaf_0
  - zlib=1.2.13=h166bdaf_4
  - zstd=1.5.2=ha4553b6_0
  - pip:
      - dssp-wsl==0.1.1
      - numpy==1.24.2
prefix: /home/labs/schwartzlab/joeg/.conda/envs/DPAM
```

