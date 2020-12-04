add_subdirectory( tools )

list( APPEND yppm_test_data
  data/inputs/yppm_0.0.1.nc
  data/inputs/yppm_0.0.3.nc
  data/inputs/yppm_0.0.6.nc
  data/inputs/yppm_0.0.8.nc
  data/inputs/yppm_0.1.3.nc
  data/inputs/yppm_0.1.6.nc
  data/inputs/yppm_0.1.7.nc
  data/inputs/yppm_0.1.15.nc
  data/inputs/yppm_0.0.7.nc
  data/inputs/yppm_0.0.9.nc
  data/inputs/yppm_0.0.12.nc
)

list( APPEND yppm_test_input
  test_input/yppm_0.0.1.nl
  test_input/yppm_0.0.3.nl
  test_input/yppm_0.0.6.nl
  test_input/yppm_0.0.8.nl
  test_input/yppm_0.1.3.nl
  test_input/yppm_0.1.6.nl
  test_input/yppm_0.1.7.nl
  test_input/yppm_0.1.15.nl
  test_input/yppm_0.0.7.nl
  test_input/yppm_0.0.9.nl
  test_input/yppm_0.0.12.nl
)

list( APPEND yppm_test_output
  test_output/yppm_0.0.1.test
  test_output/yppm_0.0.3.test
  test_output/yppm_0.0.6.test
  test_output/yppm_0.0.8.test
  test_output/yppm_0.1.3.test
  test_output/yppm_0.1.6.test
  test_output/yppm_0.1.7.test
  test_output/yppm_0.1.15.test
  test_output/yppm_0.0.7.test
  test_output/yppm_0.0.9.test
  test_output/yppm_0.0.12.test
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

# 0.0.1 input
add_test(NAME regression_0.0.1
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.1.nl 2>&1 | tee test_output/yppm_0.0.1.log.out")
add_test(NAME compare_0.0.1
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.1.log.out test_output/yppm_0.0.1.test)

# 0.0.3 input
add_test(NAME regression_0.0.3
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.3.nl 2>&1 | tee test_output/yppm_0.0.3.log.out")
add_test(NAME compare_0.0.3
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.3.log.out test_output/yppm_0.0.3.test)

# 0.0.6 input
add_test(NAME regression_0.0.6
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.6.nl 2>&1 | tee test_output/yppm_0.0.6.log.out")
add_test(NAME compare_0.0.6
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.6.log.out test_output/yppm_0.0.6.test)

# 0.0.8 input
add_test(NAME regression_0.0.8
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.8.nl 2>&1 | tee test_output/yppm_0.0.8.log.out")
add_test(NAME compare_0.0.8
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.8.log.out test_output/yppm_0.0.8.test)

# 0.1.3 input
add_test(NAME regression_0.1.3
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.1.3.nl 2>&1 | tee test_output/yppm_0.1.3.log.out")
add_test(NAME compare_0.1.3
         COMMAND ${COMPARE_SH} test_output/yppm_0.1.3.log.out test_output/yppm_0.1.3.test)

# 0.1.6 input
add_test(NAME regression_0.1.6
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.1.6.nl 2>&1 | tee test_output/yppm_0.1.6.log.out")
add_test(NAME compare_0.1.6
         COMMAND ${COMPARE_SH} test_output/yppm_0.1.6.log.out test_output/yppm_0.1.6.test)

# 0.1.7 input
add_test(NAME regression_0.1.7
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.1.7.nl 2>&1 | tee test_output/yppm_0.1.7.log.out")
add_test(NAME compare_0.1.7
         COMMAND ${COMPARE_SH} test_output/yppm_0.1.7.log.out test_output/yppm_0.1.7.test)

# 0.1.15 input
add_test(NAME regression_0.1.15
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.1.15.nl 2>&1 | tee test_output/yppm_0.1.15.log.out")
add_test(NAME compare_0.1.15
         COMMAND ${COMPARE_SH} test_output/yppm_0.1.15.log.out test_output/yppm_0.1.15.test)

# 0.0.7 input
add_test(NAME regression_0.0.7
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.7.nl 2>&1 | tee test_output/yppm_0.0.7.log.out")
add_test(NAME compare_0.0.7
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.7.log.out test_output/yppm_0.0.7.test)

# 0.0.9 input
add_test(NAME regression_0.0.9
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.9.nl 2>&1 | tee test_output/yppm_0.0.9.log.out")
add_test(NAME compare_0.0.9
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.9.log.out test_output/yppm_0.0.9.test)

# 0.0.12 input
add_test(NAME regression_0.0.12
         COMMAND bash -c "set -o pipefail; ../src/yppm test_input/yppm_0.0.12.nl 2>&1 | tee test_output/yppm_0.0.12.log.out")
add_test(NAME compare_0.0.12
         COMMAND ${COMPARE_SH} test_output/yppm_0.0.12.log.out test_output/yppm_0.0.12.test)
