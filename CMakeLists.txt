cmake_minimum_required(VERSION 3.5)
project(mitoseg)

find_package(OpenCV REQUIRED)
find_package(Boost CONFIG REQUIRED)
find_package(yaml-cpp REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})

file(GLOB SOURCES "src/*.cpp")
add_executable(mitoseg ${SOURCES})
target_link_libraries(mitoseg ${OpenCV_LIBS} yaml-cpp)
target_compile_options(mitoseg PRIVATE -Wall -Wextra -Wpedantic -Werror)
