function (cfg_compile_options OPT CFG)
  # set OPT if compiling with CFG
  add_compile_options ("$<$<CONFIG:${CFG}>:${OPT}>")
endfunction ()

include (CheckCXXCompilerFlag)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"   OR
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
    CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set (GCXX_LIKE TRUE)
  message (STATUS "GCXX-like compiler detected")
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set (ICPC_LIKE TRUE)
  message (STATUS "ICPC-like compiler detected")
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set (MSVC_LIKE TRUE)
  message (STATUS "MSVC-like compiler detected")
endif ()

if (DEFINED ICPC_LIKE)
  add_compile_options (-pedantic -Wall -Wextra -Wshadow -Wfloat-conversion)
  add_compile_options (-fp-model=precise) # needed for adapt-fp
  cfg_compile_options (-fno-math-errno RELEASE) # no errno checks
  cfg_compile_options (-fno-trapping-math RELEASE) # no fp exceptions
  cfg_compile_options (-fno-signed-zeros RELEASE) # -0.0 == +0.0
  cfg_compile_options (-ffinite-math-only RELEASE) # no explicit nans, etc
endif ()

if (DEFINED MSVC_LIKE)
  add_compile_options (/W3)
  add_compile_options (/openmp:llvm) # needed for openmp > 2.0
  add_compile_options (/fp:precise) # needed for adapt-fp
  cfg_compile_options (/GS- RELEASE) # disable buffer checks
  cfg_compile_options (/Ot RELEASE) # favour "fast" code
endif ()

if (DEFINED GCXX_LIKE)
  add_compile_options (-pedantic -Wall -Wextra -Wshadow -Wfloat-conversion)
  cfg_compile_options (-fno-math-errno RELEASE) # no errno checks
  cfg_compile_options (-fno-trapping-math RELEASE) # no fp exceptions
  cfg_compile_options (-fno-signed-zeros RELEASE) # -0.0 == +0.0
  cfg_compile_options (-ffinite-math-only RELEASE) # no explicit nans, etc
endif ()

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (INSTALL_LOCAL TRUE)
  message (STATUS "Installing locally")
else ()
  message (STATUS "Installing to ${CMAKE_INSTALL_PREFIX}")
endif ()

# try to find netcdf support

find_library (NETCDF_LIBRARY NAMES netcdf)

if (NETCDF_LIBRARY)
  message (STATUS "NetCDF library found")
  message (STATUS "NetCDF inc. lib: ${NETCDF_LIBRARY}")
else ()
  message (STATUS "NetCDF library not found")
endif ()

# try to find openmp support

find_package (OpenMP)

if (OpenMP_CXX_FOUND AND (OpenMP_CXX_VERSION LESS 3))
  set (OpenMP_CXX_FOUND FALSE)
  message (STATUS "Insufficient OpenMP support")
endif ()

if (OpenMP_CXX_FOUND)
  message (STATUS "OpenMP library found")
  message (STATUS "OpenMP inc. lib: ${OpenMP_CXX_LIB_NAMES}")
else ()
  message (STATUS "OpenMP library not found")
endif ()


add_executable (jigsaw-cmd jigsaw.cpp)
target_compile_definitions (jigsaw-cmd PRIVATE CMD_JIGSAW)

if (NETCDF_LIBRARY)
  target_compile_definitions (jigsaw-cmd PRIVATE USE_NETCDF)
  target_link_libraries (jigsaw-cmd PRIVATE ${NETCDF_LIBRARY})
  set_target_properties (jigsaw-cmd PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

if (OpenMP_CXX_FOUND)
  target_compile_definitions (jigsaw-cmd PRIVATE USE_OPENMP)
  target_link_libraries (jigsaw-cmd PRIVATE OpenMP::OpenMP_CXX)
  set_target_properties (jigsaw-cmd PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

set_target_properties (jigsaw-cmd PROPERTIES OUTPUT_NAME jigsaw)
if (DEFINED INSTALL_LOCAL)
  install (TARGETS jigsaw-cmd DESTINATION "${PROJECT_SOURCE_DIR}/bin")
else ()
  install (TARGETS jigsaw-cmd DESTINATION bin)
endif ()


add_executable (tripod-cmd jigsaw.cpp)
target_compile_definitions (tripod-cmd PRIVATE CMD_TRIPOD)

if (NETCDF_LIBRARY)
  target_compile_definitions (tripod-cmd PRIVATE USE_NETCDF)
  target_link_libraries (tripod-cmd PRIVATE ${NETCDF_LIBRARY})
  set_target_properties (tripod-cmd PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

if (OpenMP_CXX_FOUND)
  target_compile_definitions (tripod-cmd PRIVATE USE_OPENMP)
  target_link_libraries (tripod-cmd PRIVATE OpenMP::OpenMP_CXX)
  set_target_properties (tripod-cmd PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

set_target_properties (tripod-cmd PROPERTIES OUTPUT_NAME tripod)
if (DEFINED INSTALL_LOCAL)
  install (TARGETS tripod-cmd DESTINATION "${PROJECT_SOURCE_DIR}/bin")
else ()
  install (TARGETS tripod-cmd DESTINATION bin)
endif ()


add_executable (marche-cmd jigsaw.cpp)
target_compile_definitions (marche-cmd PRIVATE CMD_MARCHE)

if (NETCDF_LIBRARY)
  target_compile_definitions (marche-cmd PRIVATE USE_NETCDF)
  target_link_libraries (marche-cmd PRIVATE ${NETCDF_LIBRARY})
  set_target_properties (marche-cmd PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

if (OpenMP_CXX_FOUND)
  target_compile_definitions (marche-cmd PRIVATE USE_OPENMP)
  target_link_libraries (marche-cmd PRIVATE OpenMP::OpenMP_CXX)
  set_target_properties (marche-cmd PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

set_target_properties (marche-cmd PROPERTIES OUTPUT_NAME marche)
if (DEFINED INSTALL_LOCAL)
  install (TARGETS marche-cmd DESTINATION "${PROJECT_SOURCE_DIR}/bin")
else ()
  install (TARGETS marche-cmd DESTINATION bin)
endif ()


add_library (jigsaw-lib SHARED jigsaw.cpp)
target_compile_definitions (jigsaw-lib PRIVATE LIB_JIGSAW)

if (NETCDF_LIBRARY)
  target_compile_definitions (jigsaw-lib PRIVATE USE_NETCDF)
  target_link_libraries (jigsaw-lib PRIVATE ${NETCDF_LIBRARY})
  set_target_properties (jigsaw-lib PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

if (OpenMP_CXX_FOUND)
  target_compile_definitions (jigsaw-lib PRIVATE USE_OPENMP)
  target_link_libraries (jigsaw-lib PRIVATE OpenMP::OpenMP_CXX)
  set_target_properties (jigsaw-lib PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

set_target_properties (jigsaw-lib PROPERTIES OUTPUT_NAME jigsaw)
if (DEFINED INSTALL_LOCAL)
  install (TARGETS jigsaw-lib DESTINATION "${PROJECT_SOURCE_DIR}/lib")
else ()
  install (DIRECTORY inc/ DESTINATION include FILES_MATCHING PATTERN "*.h")
  install (TARGETS jigsaw-lib DESTINATION lib)
endif ()
