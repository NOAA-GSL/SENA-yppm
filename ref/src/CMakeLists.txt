list( APPEND yppm_src_files
netCDFModule.f90
print_affinity.c
runnable.c
yppm_core.F90
yppm_driver.f90
)

add_executable( ${PROJECT_NAME} ${yppm_src_files} )

target_include_directories( ${PROJECT_NAME} SYSTEM PUBLIC ${NetCDF_INCLUDE_DIRS} )

target_link_libraries(${PROJECT_NAME} PUBLIC NetCDF::NetCDF_Fortran)
if(OpenMP_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_C)
endif()

install(TARGETS yppm
        DESTINATION ${PROJECT_SOURCE_DIR}/exe)
