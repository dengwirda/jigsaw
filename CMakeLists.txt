cmake_minimum_required (VERSION 3.0.2)
project (JIGSAW)

include_directories ("${PROJECT_SOURCE_DIR}/inc")

add_subdirectory (src)

add_subdirectory (uni)

install (TARGETS jigsaw-bin tripod DESTINATION bin)
install (DIRECTORY inc/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install (TARGETS jigsaw DESTINATION lib)


