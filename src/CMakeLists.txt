cmake_minimum_required(VERSION 3.12)

cmake_policy(SET CMP0072 NEW) # to use GLVND

find_package(glfw3 3.3 REQUIRED)
find_package(OpenGL REQUIRED)

include_directories(
    "${CMAKE_SOURCE_DIR}/include"
    $ENV{MUJOCO_ROOT}/include
)

link_directories(
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    $ENV{MUJOCO_ROOT}/build/lib
)

file(GLOB SRCS *.cpp *.h)
add_executable(${PROJECT_NAME} ${SRCS})

target_include_directories(${PROJECT_NAME} PUBLIC ${OPENGL_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC mujoco)
target_link_libraries(${PROJECT_NAME} PUBLIC glfw Xrandr X11 Xi Xinerama ${OPENGL_gl_LIBRARY} ${CMAKE_DL_LIBS})
target_compile_options(${PROJECT_NAME} PUBLIC -Wall -Werror)