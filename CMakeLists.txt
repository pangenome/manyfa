cmake_minimum_required(VERSION 3.10)
project(ThreadSafeFaidx VERSION 0.1.0 LANGUAGES CXX)

# Build options
option(BUILD_STATIC "Build static binary" OFF)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set static compilation flags if requested
if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(CMAKE_POSITION_INDEPENDENT_CODE OFF)
else ()
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

# Create proper directory structure
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Find required packages
find_package(Threads REQUIRED)

# Find HTSlib - manual path option
option(HTSLIB_ROOT "Path to HTSlib installation" "")
if(HTSLIB_ROOT)
    include_directories(${HTSLIB_ROOT}/include)
    link_directories(${HTSLIB_ROOT}/lib)
    set(HTSLIB_LIBRARIES hts)
else()
    # Try pkg-config for system installations
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(HTSLIB QUIET htslib)
        if(HTSLIB_FOUND)
            include_directories(${HTSLIB_INCLUDE_DIRS})
            link_directories(${HTSLIB_LIBRARY_DIRS})
        else()
            # Fallback to assuming htslib is installed in a standard location
            set(HTSLIB_LIBRARIES hts)
        endif()
    else()
        # Fallback to assuming htslib is installed in a standard location
        set(HTSLIB_LIBRARIES hts)
    endif()
endif()

# Find required dependencies
find_package(ZLIB)
if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    add_definitions(-DHAVE_ZLIB)
endif()

# Find libdeflate (required for static builds with htslib)
find_library(LIBDEFLATE_LIBRARY NAMES deflate libdeflate)
if(NOT LIBDEFLATE_LIBRARY)
    message(STATUS "libdeflate not found in standard location, will try to build statically")
    if(BUILD_STATIC)
        message(WARNING "Static build enabled but libdeflate not found - static linking may fail")
    endif()
endif()

# Find additional libraries needed for static linking
find_library(LIBLZMA_LIBRARY NAMES lzma liblzma)
find_library(LIBBZ2_LIBRARY NAMES bz2 libbz2)

# Add executables
add_executable(fasta_reader src/fasta_reader.cpp)
add_executable(manyfa src/manyfa.cpp)

# Set up link libraries list based on static/dynamic build
set(COMMON_LIBRARIES
    Threads::Threads
    ${HTSLIB_LIBRARIES}
    ${ZLIB_LIBRARIES}
)

# Add additional libraries for static linking
if(BUILD_STATIC)
    list(APPEND COMMON_LIBRARIES
        ${LIBDEFLATE_LIBRARY}
        ${LIBLZMA_LIBRARY}
        ${LIBBZ2_LIBRARY}
        dl
        m
    )
endif()

# Link libraries
target_link_libraries(fasta_reader ${COMMON_LIBRARIES})
target_link_libraries(manyfa ${COMMON_LIBRARIES})

# Installation targets
install(FILES include/thread_safe_faidx.hpp DESTINATION include)
install(TARGETS fasta_reader manyfa DESTINATION bin)

# Option for building tests
option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
