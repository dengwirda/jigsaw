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

add_executable (test2d_a test2d_a.c)
target_link_libraries (test2d_a ${LIBJIGSAW})
set_target_properties(test2d_a PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2d_a DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test3d_a test3d_a.c)
target_link_libraries (test3d_a ${LIBJIGSAW})
set_target_properties(test3d_a PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test3d_a DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2s_a test2s_a.c)
target_link_libraries (test2s_a ${LIBJIGSAW})
set_target_properties(test2s_a PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2s_a DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2d_b test2d_b.c)
target_link_libraries (test2d_b ${LIBJIGSAW})
set_target_properties(test2d_b PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2d_b DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test3d_b test3d_b.c)
target_link_libraries (test3d_b ${LIBJIGSAW})
set_target_properties(test3d_b PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test3d_b DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2s_b test2s_b.c)
target_link_libraries (test2s_b ${LIBJIGSAW})
set_target_properties(test2s_b PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2s_b DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2d_c test2d_c.c)
target_link_libraries (test2d_c ${LIBJIGSAW})
set_target_properties(test2d_c PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2d_c DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test3d_c test3d_c.c)
target_link_libraries (test3d_c ${LIBJIGSAW})
set_target_properties(test3d_c PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test3d_c DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2s_c test2s_c.c)
target_link_libraries (test2s_c ${LIBJIGSAW})
set_target_properties(test2s_c PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2s_c DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2d_d test2d_d.c)
target_link_libraries (test2d_d ${LIBJIGSAW})
set_target_properties(test2d_d PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2d_d DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test3d_d test3d_d.c)
target_link_libraries (test3d_d ${LIBJIGSAW})
set_target_properties(test3d_d PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test3d_d DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2d_e test2d_e.c)
target_link_libraries (test2d_e ${LIBJIGSAW})
set_target_properties(test2d_e PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2d_e DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test3d_e test3d_e.c)
target_link_libraries (test3d_e ${LIBJIGSAW})
set_target_properties(test3d_e PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test3d_e DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2s_e test2s_e.c)
target_link_libraries (test2s_e ${LIBJIGSAW})
set_target_properties(test2s_e PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2s_e DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2d_f test2d_f.c)
target_link_libraries (test2d_f ${LIBJIGSAW})
set_target_properties(test2d_f PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2d_f DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test3d_f test3d_f.c)
target_link_libraries (test3d_f ${LIBJIGSAW})
set_target_properties(test3d_f PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test3d_f DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2s_f test2s_f.c)
target_link_libraries (test2s_f ${LIBJIGSAW})
set_target_properties(test2s_f PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2s_f DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2d_g test2d_g.c)
target_link_libraries (test2d_g ${LIBJIGSAW})
set_target_properties(test2d_g PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2d_g DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test3d_g test3d_g.c)
target_link_libraries (test3d_g ${LIBJIGSAW})
set_target_properties(test3d_g PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test3d_g DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2s_g test2s_g.c)
target_link_libraries (test2s_g ${LIBJIGSAW})
set_target_properties(test2s_g PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2s_g DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2d_h test2d_h.c)
target_link_libraries (test2d_h ${LIBJIGSAW})
set_target_properties(test2d_h PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2d_h DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2s_h test2s_h.c)
target_link_libraries (test2s_h ${LIBJIGSAW})
set_target_properties(test2s_h PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2s_h DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2d_i test2d_i.c)
target_link_libraries (test2d_i ${LIBJIGSAW})
set_target_properties(test2d_i PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2d_i DESTINATION "${PROJECT_SOURCE_DIR}")

add_executable (test2s_i test2s_i.c)
target_link_libraries (test2s_i ${LIBJIGSAW})
set_target_properties(test2s_i PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS test2s_i DESTINATION "${PROJECT_SOURCE_DIR}")
