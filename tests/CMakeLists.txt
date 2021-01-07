cmake_minimum_required(VERSION 2.8)
enable_language(CXX)
enable_language(C)

option(WITH_TCMALLOC	"Use tcmalloc from gperftools" OFF)
option(WITH_JEMALLOC	"Use jemalloc from FaceBook" OFF)

set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE "Debug")

add_definitions( "-g -O3 -Wpedantic -Wall -pthread -std=c++11" )
#add_definitions( "-DWITH_ADD_2007bl" )

include(FindPkgConfig)

enable_testing()

find_package(Boost 1.53)

if(Boost_FOUND)
	message(STATUS "Use Boost version ${Boost_VERSION}")
	include_directories(${Boost_INCLUDE_DIRS})
endif()

if (WITH_TCMALLOC)
	pkg_check_modules (LIBTCMALLOC libtcmalloc_minimal)
endif()

if (WITH_JEMALLOC AND NOT CMAKE_CROSSCOMPILING)
	pkg_check_modules (JEMALLOC jemalloc)
endif()

project(eccTest)

include_directories(.)
include_directories(../include)
#pkg_check_modules (LIBGTEST REQUIRED gtest)
pkg_check_modules (BENCH benchmark)
#find_package(GTest REQUIRED)

add_executable(ecc.test
		test.cc
		testData.hpp
        )
target_link_libraries(ecc.test gtest pthread)

if (BENCH_FOUND)
		message(STATUS "Build benchmark for vli ECC")
add_executable(ecc.bench
        bench.cc
		testData.hpp
        )
if (LIBTCMALLOC_FOUND)
	message (STATUS "Using tcmalloc...")
	target_link_libraries(ecc.bench benchmark pthread "${LIBTCMALLOC_LDFLAGS}")
elseif(JEMALLOC_FOUND)
	message (STATUS "Using jemalloc...")
	target_link_libraries(ecc.bench benchmark pthread "${JEMALLOC_LDFLAGS}")
else()
	target_link_libraries(ecc.bench benchmark pthread)
endif()
endif()

add_test(test ./ecc.test)
SET_TESTS_PROPERTIES(test PROPERTIES ENVIRONMENT "CTEST_OUTPUT_ON_FAILURE = 1")