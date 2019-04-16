set (CMAKE_CXX_STANDARD 11)
set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

add_executable (jigsaw-cmd jigsaw.cpp)
target_compile_definitions (jigsaw-cmd PRIVATE __cmd_jigsaw)
set_target_properties (jigsaw-cmd PROPERTIES OUTPUT_NAME jigsaw)
install (TARGETS jigsaw-cmd DESTINATION "${PROJECT_SOURCE_DIR}/bin")

add_executable (tripod-cmd jigsaw.cpp)
target_compile_definitions (tripod-cmd PRIVATE __cmd_tripod)
set_target_properties (tripod-cmd PROPERTIES OUTPUT_NAME tripod)
install (TARGETS tripod-cmd DESTINATION "${PROJECT_SOURCE_DIR}/bin")

add_library (jigsaw-lib SHARED jigsaw.cpp)
target_compile_definitions (jigsaw-lib PRIVATE __lib_jigsaw)
set_target_properties (jigsaw-lib PROPERTIES OUTPUT_NAME jigsaw)
install (TARGETS jigsaw-lib DESTINATION "${PROJECT_SOURCE_DIR}/lib")
