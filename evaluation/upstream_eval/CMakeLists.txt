cmake_minimum_required(VERSION 3.10)
project(score_matched_vcf)

set(CMAKE_CXX_STANDARD 17)

find_package(ZLIB REQUIRED)
find_package(CURL REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(Threads REQUIRED)

include_directories(
    $ENV{CONDA_PREFIX}/include
)

add_executable(score_matched_vcf score_matched_vcf.cpp)
target_link_libraries(score_matched_vcf 
    hts
    CURL::libcurl
    BZip2::BZip2
    LibLZMA::LibLZMA
    ZLIB::ZLIB
    Threads::Threads
)
