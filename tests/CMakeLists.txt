cmake_minimum_required(VERSION 2.8)
enable_language(CXX)
enable_language(C)

option(WITH_TCMALLOC	"Use tcmalloc from gperftools" OFF)
option(WITH_JEMALLOC	"Use jemalloc from FaceBook" OFF)
option(WITH_ENTROPY		"Use getentropy() from glibc 2.25 or above" ON)

set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE "Debug")

#add_definitions( "-Ofast -march=native -pthread -std=c++17 -funroll-loops" )
add_definitions( "-Ofast -pthread -std=c++17 -funroll-loops" )
add_definitions( "-Wpedantic -Wall -Wuninitialized" )
#add_definitions( "-DWITH_ADD_2007bl" )
#add_definitions( "-DWITH_DBL_2004hmv" )
#add_definitions( "-DWITH_MONT_D" )
add_definitions( "-DWITH_ASM" )
add_definitions( "-DNO_U64ZERO" )
# WITH_GETENTROPY only for glibc 2.25 or above (CentOS8...)

if(WITH_ENTROPY)
add_definitions( "-DWITH_GETENTROPY" )
endif()

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
