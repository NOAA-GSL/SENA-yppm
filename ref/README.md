# The YPPM reference kernel

NOTE: If you are reading this with a plain text editor, please note that this document is
formatted with Markdown syntax elements.  See https://www.markdownguide.org/cheat-sheet/
for more information.

This is the reference implementation of the `yppm` kernel extracted from the FV3 model.
This kernel does not contain any OpenMP threading. Current analysis is limited to serial, 
single threaded performance.

## Dependencies
The following packages are required for building and running this kernel:

* C compiler
* Fortran compiler
* [netcdf-c](https://www.unidata.ucar.edu/downloads/netcdf/)
* [netcdf-fortran](https://www.unidata.ucar.edu/downloads/netcdf/)
* [cmake](https://cmake.org/download/) (version >= 3.10)
* git
* [git-lfs](https://git-lfs.github.com/)

## Prerequisites
This code requires git-lfs. Before cloning the repository, verify that
git-lfs is installed, by issuing the following command. This only needs to be done once per user per machine.

```bash
$ git lfs install
```

If the above gives an error you (or your systems administrator) may need to install git-lfs.

Some systems that use modules to manage software provide git with git-lfs support via a
module (e.g. `module load git`).  If you are using a system that uses modules, use
`module avail` to look for alternative versions of git that may have git-lfs support.

Make sure the files in `data/inputs` are NetCDF data files (not text) before proceeding to
the build step. A simple way to do that is with the `file` command as shown below:

```
$ file data/inputs/*
inputs/yppm_0.0.1.nc:  NetCDF Data Format data
inputs/yppm_0.0.12.nc: NetCDF Data Format data
inputs/yppm_0.0.3.nc:  NetCDF Data Format data
inputs/yppm_0.0.6.nc:  NetCDF Data Format data
inputs/yppm_0.0.7.nc:  NetCDF Data Format data
inputs/yppm_0.0.8.nc:  NetCDF Data Format data
inputs/yppm_0.0.9.nc:  NetCDF Data Format data
inputs/yppm_0.1.15.nc: NetCDF Data Format data
inputs/yppm_0.1.3.nc:  NetCDF Data Format data
inputs/yppm_0.1.6.nc:  NetCDF Data Format data
inputs/yppm_0.1.7.nc:  NetCDF Data Format data
```

**NOTE**: If you cloned the repository with a version of git without git-lfs installed,
or before you ran `git lfs install` you
must run the following command (with a version of git that does support git-lfs) from the base
of the repository to fetch the input data before proceeding to the build steps. 

```bash
$ git lfs pull
```

Alternatively, you can reclone the repository with git-lfs installed.

## Building the kernel

This kernel uses an out-of-source cmake build, meaning that the build must be done in 
directory that is not in the source tree.

### Basic build procedure (from the directory containing this file)

The basic build steps are as follows (**NOTE**: Make sure not to omit the two dots at the end
of the `cmake` step.):

```bash
$ rm -rf build ; mkdir build ; cd build
$ export CC=<name of C compiler>
$ export FC=<name of fortran compiler> 
$ cmake -DCMAKE_BUILD_TYPE=<debug | release> ..
$ make VERBOSE=1
```

On many systems, the above will suffice. However, some systems will require you to help cmake
find dependencies, particularly if software depencencies are not installed in standard locations.
See below for more information.

### Machines that use modules to manage software

Most HPC systems use modules to manage software.  Make sure you have loaded the versions of
the compiler and software you want to use before running the build steps above.  This will allow build
dependencies to be found properly.  For example:

```bash
$ module load intel netcdf cmake
```

### Machines that do not use modules to manage software

If compilers and/or NetCDF is not installed in a standard location where cmake can find it, you
may need to add their installation paths to the `CMAKE_PREFIX_PATH` before running the steps
above. For example:

```bash
$ export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/path/to/netcdf:/path/to/netcdf-fortran
```

### Building on a Mac

By default, gcc points to the clang compiler on Mac.  To use the GNU compiler on Mac, depending
on how the GNU compiler was installed, you may need to specify the C compiler name as gcc-$version.
For example:

```bash
$ export CC=gcc-10
```

## Testing the kernel

First, set the OpenMP variable. Although `yppm` is currently single-threaded, it may be threaded in the future, and setting `OMP_NUM_THREADS=1` avoids extraneous output of thread affinity information. 

```bash
$ export OMP_NUM_THREADS=1
```

Then, to run the test suite (from the build directory):

```bash
$ ctest
```

To run a specific test (for example):

```bash
$ ctest -R regression_0.0.3
```

To run a specific test with full output to get more information about a failure (for example):

```bash
$ ctest -VV -R regression_0.0.3
```

## Build and test script

For convenience, a build script is provided that builds the code and runs the test suite.

**(NOTE: This script is written for machines that use modules and it may need to be modified,
depending on how modules are set up on your system)**

```bash
$ ./build.sh <gcc | intel> <debug | release>
```

## Installation and running

To (optionally) install the built executable into exe/

```bash
$ make install
```

To run the installed executable (for example):

```bash
$ export OMP_NUM_THREADS=1
$ exe/yppm ../test/test_input/yppm_0.0.1.nl
```

## NOTES:

1. The test suite does not measure performance, but reports how long each test takes to run.
2. Detailed performance timings are printed in the stdout when the kernel runs.
3. To view kernel output, either run ctest in verbose mode, or run the kernel manually.

## Here is a list of the files and what they contain.

- `cmake/` contains compiler flag settings and cmake helper scripts
- `src/` contains the kernel source code
- `test/` contains the tests, test input, and test output
- `test/data/outputs` is where test output data is written
- `test/test_input` contains the test namelist input files
- `test/test_output` contains the test baselines
- `exe/` contains the installed executable

## Troubleshooting

1. All tests fail on my machine.

    Check to make sure git-lfs is installed and that all files in `data/inputs` are NetCDF 
    data files and are not text. Run `git lfs pull` to download NetCDF files if necessary.

2. I get `Skipping object checkout, Git LFS is not installed.` when running `git lfs pull`

    Run `git lfs install` to perform the one-time installation that git-lfs requires per user per machine.

3. I get `git: 'lfs' is not a git command.` when running `git lfs pull`

    Your version of git does not support git-lfs. Install git-lfs or load a version of git that supports it.

4. I get `git-lfs smudge -- 'data/inputs/yppm_0.0.1.nl': git-lfs: command not found` when cloning.

    Your version of git does not support git-lfs. Install git-lfs or load a version of git that supports it.

5. I get unresolved symbols when testing / running the kernel

    If you are on a machine that uses modules to manage software, you probably need to load the modules
    for your compiler and/or NetCDF **(make sure to use the same modules to build, test, and run)**.  For example:
    ```bash
    $ module load intel netcdf
    ```

    If you are on a machine that does not use modules, you probably need to add the paths of your compiler
    and/or NetCDF libraries to `LD_LIBRARY_PATH`.  For example:
    ```bash
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/netcdf-c/lib:/path/to/netcdf-fortran/lib
    ```

6. Executing the kernels, I get "forrtl: severe (174): SIGSEGV, segmentation fault occurred"
 
    This error may be due to too small of a stack size.  Try increasing the stack size with:
    ```bash
    $ limit -s 5000000 
    ```
