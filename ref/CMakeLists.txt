################################################################################
# yppm
################################################################################

cmake_minimum_required( VERSION 3.10 )

project( yppm LANGUAGES Fortran C )

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

enable_testing()

include(${PROJECT_NAME}_compiler_flags)

find_package( OpenMP COMPONENTS C Fortran )
find_package( NetCDF REQUIRED COMPONENTS C Fortran )

add_subdirectory(src)
add_subdirectory(test)
