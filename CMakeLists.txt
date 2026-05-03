cmake_minimum_required(VERSION 3.10)
project(PricerLib VERSION 1.0)

# Force C++17 for cross-platform std::filesystem support
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include current directory to find headers like "RandomGenerator/Normal.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ---------------------------------------------------------------------------
# Core Library
# ---------------------------------------------------------------------------
# We collect all .cpp files from the core directories
file(GLOB_RECURSE CORE_SOURCES 
    "RandomGenerator/*.cpp"
    "SDE/*.cpp"
    "Payoffs/*.cpp"
    "Pricer/*.cpp"
    "VarRed/*.cpp"
)

# Exclude Visual Studio specific files and template/main files if any
list(FILTER CORE_SOURCES EXCLUDE REGEX "pch.cpp|dllmain.cpp|main.cpp")

add_library(PricerCore STATIC ${CORE_SOURCES})

# ---------------------------------------------------------------------------
# Executables (Demos & Studies)
# ---------------------------------------------------------------------------
add_executable(run_demo demo/demo.cpp)
target_link_libraries(run_demo PricerCore)

add_executable(run_study demo/study.cpp)
target_link_libraries(run_study PricerCore)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
# Create an executable for each test file in the tests/ directory
file(GLOB TEST_SOURCES "tests/test_*.cpp")

foreach(TEST_FILE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE})
    target_link_libraries(${TEST_NAME} PricerCore)
endforeach()

# ---------------------------------------------------------------------------
# Output Directory settings (optional but clean)
# ---------------------------------------------------------------------------
set_target_properties(run_demo run_study PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
