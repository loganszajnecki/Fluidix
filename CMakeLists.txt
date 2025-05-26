cmake_minimum_required(VERSION 3.15)
project(Fluidix LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Output directories for binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include paths
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/core
    ${PROJECT_SOURCE_DIR}/include/flux
    ${PROJECT_SOURCE_DIR}/include/solver
    ${PROJECT_SOURCE_DIR}/include/reconstruction
    ${PROJECT_SOURCE_DIR}/include/io
)

# Source groups for IDEs
file(GLOB_RECURSE CORE_SRC src/core/*.cpp)
file(GLOB_RECURSE FLUX_SRC src/flux/*.cpp)
file(GLOB_RECURSE SOLVER_SRC src/solver/*.cpp)
file(GLOB_RECURSE RECON_SRC src/reconstruction/*.cpp)
file(GLOB_RECURSE IO_SRC src/io/*.cpp)

# Group all together
set(SOURCES
    ${CORE_SRC}
    ${FLUX_SRC}
    ${SOLVER_SRC}
    ${RECON_SRC}
    ${IO_SRC}
    main.cpp
)

# Executable
add_executable(Fluidix ${SOURCES})
