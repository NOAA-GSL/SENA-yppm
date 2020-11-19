if( NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_definitions( -DNDEBUG )
endif( )

# Set flags depending on which compiler iis used
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  include("compiler_flags_GNU_Fortran")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  include("compiler_flags_Intel_Fortran")
endif()

# Set flags depending on which compiler iis used
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  include("compiler_flags_GNU_C")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  include("compiler_flags_Intel_C")
endif()
