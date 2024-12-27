#!/bin/bash 

echo "Unloading modules"
module unload llvm-14.0.0
module unload nvhpc-23.9
module unload cuda-11.4.4

# Check for a command-line argument
if [ "$1" == "llvm" ]; then
    echo "Loading clang modules"
    module load llvm-14.0.0
    module load cuda-11.4.4
elif [ "$1" == "nvhpc" ]; then
    echo "Loading nvhpc modules"
    module load nvhpc-23.9
else
    echo "Usage: $0 [llvm|nvhpc]"
    exit 1
fi

module load cmake-3.23.0
module load petsc-3.17.1
module load openmpi-4.0.2
