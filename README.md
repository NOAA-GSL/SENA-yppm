![CI for yppm](https://github.com/NOAA-GSL/SENA-yppm/workflows/CI%20for%20yppm/badge.svg?branch=develop)

```
This repository is a scientific product and is not official communication
of the National Oceanic and Atmospheric Administration, or the United States
Department of Commerce. All NOAA GitHub project code is provided on an ‘as
is’ basis and the user assumes responsibility for its use. Any claims against
the Department of Commerce or Department of Commerce bureaus stemming from
the use of this GitHub project will be governed by all applicable Federal
law. Any reference to specific commercial products, processes, or service
by service mark, trademark, manufacturer, or otherwise, does not constitute
or imply their endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and logo of
a DOC bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.
```

This project contains derivative work from the Geophysical Fluid Dynamics
Laboratory GitHub project for the FV3 dynamical core
(https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere). It is a stand alone
subset of the original specifically licensed for use under the
[Apache License](https://www.apache.org/licenses/LICENSE-2.0).
The application of the Apache License to this standalone work does not imply
nor shall it be construed as any change to the original source license.

# Overview

This repository contains a stand-alone kernel for the `yppm` subroutine,
extracted from the FV3 model. The goal of this repository is to facilitate
direct comparisons of different implementations of `yppm`. Alternative
implementations may use different programming languages, parallel programming
models, or numerical modeling libraries. They may also target different
hardware such as GPUs or ARM processors. In numerical modeling, metrics
for comparison traditionally include performance and portability. However,
this work encourages comparison of developer productivity and design metrics
such as ease of use, extensibility, and maintainability.

A Fortran 90 implementation of `yppm` is provided as a baseline reference
to which other implementations should be compared. Additional implementations
will be added as they are developed.

# Contents

This repository is organized as follows.

### `data/`

Contains the reference input and output data all implementations must use
for testing and validation.

### `ref/`

Contains the source tree for the reference implementation. See the README
in that directory for details about how to run and test it.

### `foo/`

Alternative implementations will be stored in a new top-level subdirectory
along side the reference implementation. In this documentation, `foo/` is
simply a place holder for names given to future implementations.

# Contributing

Please see the [Contributing Guide](CONTRIBUTING.md) for information about contributing to this repository.

