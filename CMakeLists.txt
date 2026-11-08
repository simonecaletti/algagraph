cmake_minimum_required(VERSION 3.14)
project(ALGAgraph VERSION 2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(algacore STATIC src/linalg.cpp src/alga.cpp src/plot.cpp)
target_include_directories(algacore PUBLIC src)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(algacore PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(algagraph src/main.cpp)
target_link_libraries(algagraph PRIVATE algacore)

include(CTest)
if(BUILD_TESTING)
  add_executable(test_alga tests/test_alga.cpp)
  target_link_libraries(test_alga PRIVATE algacore)
  add_test(NAME unit COMMAND test_alga)
endif()

install(TARGETS algagraph RUNTIME DESTINATION bin)
