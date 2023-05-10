cmake_minimum_required(VERSION 3.8)

include_directories(
    "${CMAKE_SOURCE_DIR}/include"
    $ENV{MUJOCO_ROOT}/include
)

link_directories(
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    $ENV{MUJOCO_ROOT}/lib
)

file(GLOB SRCS *.cpp *.h)
add_executable(${PROJECT_NAME} ${SRCS})
target_link_libraries(
    ${PROJECT_NAME}
    mujoco
)