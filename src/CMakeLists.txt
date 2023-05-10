cmake_minimum_required(VERSION 3.8)

find_package(glfw3 3.3 REQUIRED)
find_package(OpenGL REQUIRED)

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
    PUBLIC
    mujoco glfw3
)
target_include_directories(${PROJECT_NAME} PUBLIC ${OPENGL_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC glfw ${OPENGL_gl_LIBRARY} ${CMAKE_DL_LIBS})