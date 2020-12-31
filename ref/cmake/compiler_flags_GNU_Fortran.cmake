####################################################################
# COMMON FLAGS
####################################################################
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -ffloat-store -fcray-pointer -fno-unsafe-math-optimizations -frounding-math -mno-fused-madd")

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native -funroll-all-loops -finline-functions")

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O1 -g -fcheck=bounds -ffpe-trap=invalid,zero,overflow" )
