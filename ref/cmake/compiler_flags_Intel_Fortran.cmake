####################################################################
# COMMON FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback -i4 -real-size 32 -fpp -assume byterecl -safe-cray-ptr -ftz -fp-model source" )

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-fno-alias -auto -align array64byte -xHOST -qno-opt-dynamic-align -O2 -qoverride-limits -qopt-prefetch=3" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0 -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -init=snan,array -warn unused" )
