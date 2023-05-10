cmake_minimum_required(VERSION 3.8)

include_directories(
    "${CMAKE_SOURCE_DIR}/include"
)

#link_directories(
#    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
#)

file(GLOB SRCS *.cpp *.h)
add_executable(${PROJECT_NAME} ${SRCS})
#target_link_libraries(
#    ${PROJECT_NAME}
#    gtest pthread
#)