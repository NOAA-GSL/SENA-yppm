#!/bin/bash

set -e

# Check usage
if [ $# -lt 2 -o $# -gt 3 ]; then
  echo "Usage: build.sh gcc|intel debug|release [off|manual|auto]"
  exit 1
fi

# Get arguments
COMPILER=$1
BUILD_TYPE=$2
GPTL="${3:-off}"

# Start with a clean module environment
module purge

# Handle compiler settings
case ${COMPILER} in
  "gcc"|"gnu" )
    export CC=gcc
    export FC=gfortran
    module load gcc 2>/dev/null || module load gnu
  ;;
  "intel" )
    export CC=icc
    export FC=ifort
    module load intel
  ;;
  *) 
    echo "Unsupported compiler: ${COMPILER}"
    exit 1
  ;;
esac

# Handle build type
case ${BUILD_TYPE} in
  "debug"|"release" )
    echo "Building with ${COMPILER} in ${BUILD_TYPE} mode"
  ;;
  *)
    echo "Unsupported build type: ${BUILD_TYPE}"
    exit 1
  ;;
esac

# Handle profile flag
GPTL_CMAKE_OPTIONS=""
case ${GPTL} in
  "off" )
    echo "Building without GPTL enabled"
  ;;
  "manual" )
    echo "Building with GPTL enabled"
    GPTL_CMAKE_OPTIONS="-DENABLE_GPTL=1"
    module load gptl
  ;;
  "auto" )
    echo "Building with GPTL auto-profiling enabled"
    GPTL_CMAKE_OPTIONS="-DENABLE_GPTL=1 -DENABLE_AUTOPROFILING=1"
    module load gptl
  ;;
  *)
    echo "Unsupportd GPTL profile mode: ${GPTL}"
    exit 1
  ;;
esac

# Print out modules in use
module load netcdf cmake
module list

# Create build dir
rm -rf build
mkdir build
cd build

# Build
cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ${GPTL_CMAKE_OPTIONS} ..
make -j4 VERBOSE=1

# Run the test suite
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close
ulimit -s 5000000

ctest
