cmake_minimum_required (VERSION 3.9.4)
project (JIGSAW_UNIT_TEST)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (INSTALL_LOCAL TRUE)
  message ( ">> INSTALLING LOCALLY" )
else ()
  message ( ">> INSTALLING TO ${CMAKE_INSTALL_PREFIX}" )
endif ()

if (DEFINED INSTALL_LOCAL)
  find_library (LIBJIGSAW jigsaw HINTS "${PROJECT_SOURCE_DIR}/../lib")
else ()
  find_library (LIBJIGSAW jigsaw HINTS "${CMAKE_INSTALL_PREFIX}/lib")
endif ()

message (">> LIBJIGSAW LOCATION: ${LIBJIGSAW}")

add_executable (test_1 test_1.c)
target_link_libraries (test_1 ${LIBJIGSAW})
set_target_properties(test_1 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test_1 DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test_2 test_2.c)
target_link_libraries (test_2 ${LIBJIGSAW})
set_target_properties(test_2 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test_2 DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test_3 test_3.c)
target_link_libraries (test_3 ${LIBJIGSAW})
set_target_properties(test_3 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test_3 DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test_4 test_4.c)
target_link_libraries (test_4 ${LIBJIGSAW})
set_target_properties(test_4 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test_4 DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test_5 test_5.c)
target_link_libraries (test_5 ${LIBJIGSAW})
set_target_properties(test_5 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test_5 DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test_6 test_6.c)
target_link_libraries (test_6 ${LIBJIGSAW})
set_target_properties(test_6 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test_6 DESTINATION "${PROJECT_SOURCE_DIR}")
