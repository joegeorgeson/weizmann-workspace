This is for a GROMACS implentation of [relax](https://git.embl.de/grp-kosinski/relax/-/tree/main)

I created a conda environment and the below modules

ml Python/3.9.6-GCCcore-11.2.0 gromacs/2021.2 gromacs/2021.2 TensorFlow/2.7.1-foss-2021b-CUDA-11.4.1


```
conda create --name relax python=3.9.6
git clone https://git.embl.de/grp-kosinski/relax.git
conda activate relax
ml restore relax

python setup.py build #didn't install so use absolute paths

cd /home/labs/schwartzlab/joeg/github/relax
#this worked
bsub -R rusage[mem=1000] -gpu num=5:j_exclusive=yes -q gpu-long -C 0 -e err -o out python scripts/relax_gromacs.py --mode gpu --gpu_count 5 test/data/unrelaxed_AlphaFold_model.pdb test/data/unrelaxed_AlphaFold_model.min.pdb
```


For backup;

```
(relax) [joeg@epitrans relax]$ ml

Currently Loaded Modules:
  1) GCCcore/11.2.0                  11) libffi/3.4.2-GCCcore-11.2.0       21) libevent/2.1.12-GCCcore-11.2.0  31) CUDA/11.4.1                                 41) cURL/7.78.0-GCCcore-11.2.0              51) protobuf/3.17.3-GCCcore-11.2.0
  2) zlib/1.2.11-GCCcore-11.2.0      12) OpenSSL/1.1                       22) UCX/1.11.2-GCCcore-11.2.0       32) cuDNN/8.2.2.26-CUDA-11.4.1                  42) double-conversion/3.1.5-GCCcore-11.2.0  52) protobuf-python/3.17.3-GCCcore-11.2.0
  3) binutils/2.37-GCCcore-11.2.0    13) Python/3.9.6-GCCcore-11.2.0       23) PMIx/4.1.0-GCCcore-11.2.0       33) GDRCopy/2.3-GCCcore-11.2.0                  43) flatbuffers/2.0.0-GCCcore-11.2.0        53) flatbuffers-python/2.0-GCCcore-11.2.0
  4) bzip2/1.0.8-GCCcore-11.2.0      14) gcc/8.3.0                         24) OpenMPI/4.1.1-GCC-11.2.0        34) UCX-CUDA/1.11.2-GCCcore-11.2.0-CUDA-11.4.1  44) giflib/5.2.1-GCCcore-11.2.0             54) libpng/1.6.37-GCCcore-11.2.0
  5) ncurses/6.2-GCCcore-11.2.0      15) gromacs/2021.2                    25) OpenBLAS/0.3.18-GCC-11.2.0      35) NCCL/2.10.3-GCCcore-11.2.0-CUDA-11.4.1      45) ICU/69.1-GCCcore-11.2.0                 55) snappy/1.1.9-GCCcore-11.2.0
  6) libreadline/8.1-GCCcore-11.2.0  16) GCC/11.2.0                        26) FlexiBLAS/3.0.4-GCC-11.2.0      36) pybind11/2.7.1-GCCcore-11.2.0               46) JsonCpp/1.9.4-GCCcore-11.2.0            56) networkx/2.6.3-foss-2021b
  7) Tcl/8.6.11-GCCcore-11.2.0       17) numactl/2.0.14-GCCcore-11.2.0     27) gompi/2021b                     37) SciPy-bundle/2021.10-foss-2021b             47) NASM/2.15.05-GCCcore-11.2.0             57) TensorFlow/2.7.1-foss-2021b-CUDA-11.4.1
  8) SQLite/3.36-GCCcore-11.2.0      18) libxml2/2.9.10-GCCcore-11.2.0     28) FFTW/3.3.10-gompi-2021b         38) Szip/2.1.1-GCCcore-11.2.0                   48) libjpeg-turbo/2.0.6-GCCcore-11.2.0
  9) XZ/5.2.5-GCCcore-11.2.0         19) libpciaccess/0.16-GCCcore-11.2.0  29) ScaLAPACK/2.1.0-gompi-2021b-fb  39) HDF5/1.12.1-gompi-2021b                     49) LMDB/0.9.29-GCCcore-11.2.0
 10) GMP/6.2.1-GCCcore-11.2.0        20) hwloc/2.5.0-GCCcore-11.2.0        30) foss/2021b                      40) h5py/3.6.0-foss-2021b                       50) nsync/1.24.0-GCCcore-11.2.0

```
