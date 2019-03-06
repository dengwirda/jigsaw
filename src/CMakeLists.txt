add_executable (jigsaw-bin jigsaw.cpp)
target_compile_definitions (jigsaw-bin PRIVATE NDEBUG __cmd_jigsaw)
set_target_properties (jigsaw-bin PROPERTIES OUTPUT_NAME jigsaw)

add_executable (tripod jigsaw.cpp)
target_compile_definitions (tripod PRIVATE NDEBUG __cmd_tripod)

add_library (jigsaw SHARED jigsaw.cpp)
target_compile_definitions (jigsaw PRIVATE NDEBUG __lib_jigsaw)
