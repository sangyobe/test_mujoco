cmake_minimum_required(VERSION 3.12)
project("test_${CMAKE_PROJECT_NAME}")

cmake_policy(SET CMP0074 NEW)

SET(EIGEN3_ROOT_DIR "$ENV{EIGEN3_ROOT}")
IF(NOT EIGEN3_ROOT_DIR)
    MESSAGE(FATAL_ERROR "Please point the environment variable EIGEN3_ROOT to the include directory of your Eigen3 installation.")
ENDIF()
find_package(Eigen3 3.4.0 REQUIRED)
if(NOT Eigen3_FOUND)
    MESSAGE(FATAL_ERROR "Eigen3 not found.")
endif()
MESSAGE(STATUS "Eigen3 was found at ${EIGEN3_INCLUDE_DIR}.")
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.13.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


include_directories(
    "${CMAKE_SOURCE_DIR}/include"
)
link_directories(
    "${CMAKE_SOURCE_DIR}/lib"
)

#link_directories(
#    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
#)

file(GLOB SRCS test_quadip*.cpp test_quadip*.h)
add_executable(${PROJECT_NAME} ${SRCS})
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} gtest gtest_main pthread)

file(GLOB dtMath_SRCS test_dtmath*.cpp test_dtmath*.h)
add_executable(test_dtMath ${dtMath_SRCS})
add_test(NAME test_dtMath COMMAND test_dtMath)
target_link_libraries(test_dtMath gtest gtest_main pthread)

message(STATUS "CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")

file(GLOB BLOB_SRCS test_blob*.cpp)
add_executable(test_blob ${BLOB_SRCS})
add_test(NAME test_blob COMMAND test_blob WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
target_link_libraries(test_blob gtest gtest_main pthread tinyxml2)