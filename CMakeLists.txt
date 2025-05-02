cmake_minimum_required(VERSION 3.10)
project(ThreadSafeFaidx VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

# Find zlib
find_package(ZLIB)
if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    add_definitions(-DHAVE_ZLIB)
endif()

# Add executables
add_executable(fasta_reader src/fasta_reader.cpp)
add_executable(manyfasta src/manyfasta.cpp)

# Link libraries
target_link_libraries(fasta_reader
    Threads::Threads
    ${HTSLIB_LIBRARIES}
    ${ZLIB_LIBRARIES}
)

target_link_libraries(manyfasta
    Threads::Threads
    ${HTSLIB_LIBRARIES}
    ${ZLIB_LIBRARIES}
)

# Installation targets
install(FILES include/thread_safe_faidx.hpp DESTINATION include)
install(TARGETS fasta_reader manyfasta DESTINATION bin)

# Option for building tests
option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
