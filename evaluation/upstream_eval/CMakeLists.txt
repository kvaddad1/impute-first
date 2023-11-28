cmake_minimum_required(VERSION 3.0)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

project (score_matched_vcf)

find_library(LIB_HTS hts)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra --std=c++17 -D__STDC_FORMAT_MACROS")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-g -Ofast -fstrict-aliasing -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

include_directories($ENV{CONDA_PREFIX}/include)
link_directories($ENV{CONDA_PREFIX}/lib)

message("Building in ${CMAKE_BUILD_TYPE} mode")

add_executable(score_matched_vcf score_matched_vcf.cpp)
target_link_libraries(score_matched_vcf z)
target_link_libraries(score_matched_vcf hts)
