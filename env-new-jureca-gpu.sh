#!/bin/bash
## MODULES
# module --force purge
module load Stages/2024  GCC/12.3.0
module load OpenMPI/4.1.5
module load MPI-settings/CUDA
module load UCX-settings/RC-CUDA
module load CMake/3.26.3
module load HDF5/1.14.2
module load GMP/6.2.1
module load Python/3.11.3
module load OpenBLAS/0.3.23
## COMPILE ENV
export GPU_ARCH=sm_80
export CC=mpicc
export CXX=mpicxx
export CMAKE=cmake
export CMAKE_EXTRA_FLAGS="-DCMAKE_EXPORT_COMPILE_COMMANDS=1"
export CMAKE_MAKE_FLAGS="--parallel $(nproc)"

## RUNTIME ENV
export CODE_BASE_DIR="/p/scratch/exotichadrons/chroma-distillation"
export LD_LIBRARY_PATH="$CODE_BASE_DIR/install/chroma/lib:$CODE_BASE_DIR/install/quda/lib:$CODE_BASE_DIR/install/qdpjit/lib:$CODE_BASE_DIR/install/qmp/lib:$CODE_BASE_DIR/install/llvm/lib:$LD_LIBRARY_PATH"
export QUDA_RESOURCE_PATH="$CODE_BASE_DIR/quda_resources"
export QUDA_TUNE_VERSION_CHECK=0

