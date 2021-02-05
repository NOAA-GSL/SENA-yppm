add_subdirectory( tools )

list( APPEND yppm_test_data
  data/inputs/yppm_0.1.7.nc
  data/inputs/yppm_0.0.12.nc
)

list( APPEND yppm_test_input
  test_input/yppm_0.0.12.nl
  test_input/yppm_0.0.12_10.nl
  test_input/yppm_0.0.12_20.nl
  test_input/yppm_0.0.12_30.nl
  test_input/yppm_0.1.7.nl
  test_input/yppm_0.1.7_10.nl
  test_input/yppm_0.1.7_20.nl
  test_input/yppm_0.1.7_30.nl
)

list( APPEND yppm_test_output
  test_output/yppm_0.0.12.test
  test_output/yppm_0.0.12_10.test
  test_output/yppm_0.0.12_20.test
  test_output/yppm_0.0.12_30.test
  test_output/yppm_0.1.7.test
  test_output/yppm_0.1.7_10.test
  test_output/yppm_0.1.7_20.test
  test_output/yppm_0.1.7_30.test
)

# Create data directory for kernel input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data/inputs)
foreach(FILENAME ${yppm_test_data})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/../../${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create empty data directory for model test output
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data/outputs)

# Create test_input directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_input)
foreach(FILENAME ${yppm_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_output)
foreach(FILENAME ${yppm_test_output})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

set(COMPARE_SH ${CMAKE_BINARY_DIR}/bin/yppm_compare.sh)

# Add unit tests here.
list (APPEND test_src_files
${PROJECT_SOURCE_DIR}/test/src/interpolateTest.f90
)
add_executable (interpolateTest ${test_src_files}  ${PROJECT_SOURCE_DIR}/src/interpolate.f90)

if ( GPTL_FOUND )
   target_link_libraries( interpolateTest PUBLIC GPTL::GPTL )
endif()

# 0.1.7 input
add_test(NAME regression_0.1.7
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.1.7.nl 2>&1 | tee test_output/yppm_0.1.7.log.out")
add_test(NAME compare_0.1.7
         COMMAND ${COMPARE_SH} test_output/yppm_0.1.7.log.out test_output/yppm_0.1.7.test)

# 0.1.7 input_10
add_test(NAME regression_0.1.7_10
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.1.7_10.nl 2>&1 | tee test_output/yppm_0.1.7_10.log.out")
add_test(NAME compare_0.1.7_10
         COMMAND ${COMPARE_SH} test_output/yppm_0.1.7_10.log.out test_output/yppm_0.1.7_10.test)

# 0.1.7 input_20
add_test(NAME regression_0.1.7_20
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.1.7_20.nl 2>&1 | tee test_output/yppm_0.1.7_20.log.out")
add_test(NAME compare_0.1.7_20
         COMMAND ${COMPARE_SH} test_output/yppm_0.1.7_20.log.out test_output/yppm_0.1.7_20.test)

# 0.1.7 input_30
add_test(NAME regression_0.1.7_30
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.1.7_30.nl 2>&1 | tee test_output/yppm_0.1.7_30.log.out")
add_test(NAME compare_0.1.7_30
         COMMAND ${COMPARE_SH} test_output/yppm_0.1.7_30.log.out test_output/yppm_0.1.7_30.test)
 
# 0.0.12 input
add_test(NAME regression_0.0.12
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.12.nl 2>&1 | tee test_output/yppm_0.0.12.log.out")
add_test(NAME compare_0.0.12
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.12.log.out test_output/yppm_0.0.12.test)

# 0.0.12 input_10
add_test(NAME regression_0.0.12_10
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.12_10.nl 2>&1 | tee test_output/yppm_0.0.12_10.log.out")
add_test(NAME compare_0.0.12_10
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.12_10.log.out test_output/yppm_0.0.12_10.test)

# 0.0.12 input_20
add_test(NAME regression_0.0.12_20
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.12_20.nl 2>&1 | tee test_output/yppm_0.0.12_20.log.out")
add_test(NAME compare_0.0.12_20
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.12_20.log.out test_output/yppm_0.0.12_20.test)

# 0.0.12 input_30
add_test(NAME regression_0.0.12_30
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.12_30.nl 2>&1 | tee test_output/yppm_0.0.12_30.log.out")
add_test(NAME compare_0.0.12_30
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.12_30.log.out test_output/yppm_0.0.12_30.test)

add_test(NAME InterpolateTest 
         COMMAND bash -c "set -o pipefail; ../test/interpolateTest 2>&1 | tee test_output/interpolateTest.out")

# Test install
add_test(NAME install
         COMMAND make install)
