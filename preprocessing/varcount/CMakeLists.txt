cmake_minimum_required(VERSION 3.10)

project(varcount CXX C)

set(PROJECT_URL "https://github.com/alshai/varcount")
set(PROJECT_DESCRIPTION "Outputs alignment coverage over a predefined set of variants")
set(CMAKE_CXX_STANDARD 11)

# htslib
find_library(LIB_HTS hts)
include_directories($ENV{CONDA_PREFIX}/include)
link_directories($ENV{CONDA_PREFIX}/lib)

add_executable(varcount varcount.cpp hts_util.hpp)
target_link_libraries(varcount ${LIB_HTS})

add_executable(vcf_score vcf_score.cpp hts_util.hpp)
target_link_libraries(vcf_score ${LIB_HTS})

INSTALL(TARGETS varcount DESTINATION bin)
