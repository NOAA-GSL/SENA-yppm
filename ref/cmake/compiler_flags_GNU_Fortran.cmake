####################################################################
# COMMON FLAGS
####################################################################
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace -ffloat-store -fcray-pointer -fno-unsafe-math-optimizations -frounding-math -ffp-contract=off")

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-all-loops -finline-functions")

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow" )

####################################################################
# FLAGS FOR AUTOPROFILING
####################################################################

set( Fortran_AUTOPROFILING_FLAGS        "-finstrument-functions" )
