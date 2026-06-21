cmake_minimum_required(VERSION 3.20)
project(moire_cpp_tests)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler-specific flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
endif()

# Find required packages
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)

# Include directories (mirror src/Makevars module layout)
include_directories(../src)
include_directories(../src/multivector)
include_directories(../src/observation)
include_directories(../src/inference)
include_directories(../src/mcmc)
include_directories(../src/core)
include_directories(support)
include_directories(.)
include_directories(${TBB_INCLUDE_DIRS})

execute_process(
    COMMAND Rscript -e "cat(system.file('include', package='BH'))"
    OUTPUT_VARIABLE BH_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if (BH_INCLUDE_DIR)
    include_directories(${BH_INCLUDE_DIR})
endif()

# Shared parallel support (TBB via multivector_ops.h; no OpenMP)
add_library(moire_parallel_support INTERFACE)
target_link_libraries(moire_parallel_support INTERFACE Threads::Threads TBB::tbb)

# Test framework (Catch2)
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.4.0
)
FetchContent_MakeAvailable(Catch2)

# Test executables
add_executable(multivector_tests tests/multivector_tests.cpp)
target_link_libraries(multivector_tests
    PRIVATE
    Catch2::Catch2WithMain
    moire_parallel_support
)

add_executable(distance_matrix_tests tests/distance_matrix_test.cpp)
target_link_libraries(distance_matrix_tests
    PRIVATE
    Catch2::Catch2WithMain
    moire_parallel_support
)

add_executable(observation_model_tests tests/observation_model_tests.cpp)
target_link_libraries(observation_model_tests
    PRIVATE
    Catch2::Catch2WithMain
    moire_parallel_support
)
target_sources(observation_model_tests PRIVATE
    tests/observation_model_test_stubs.cpp
)

add_executable(prob_any_missing_tests tests/prob_any_missing_tests.cpp)
target_link_libraries(prob_any_missing_tests
    PRIVATE
    Catch2::Catch2WithMain
    moire_parallel_support
)
target_sources(prob_any_missing_tests PRIVATE
    ../src/prob_any_missing.cpp
    ../src/combination_indices_generator.cpp
)

add_executable(zeta_mobius_tests tests/zeta_mobius_tests.cpp)
target_link_libraries(zeta_mobius_tests
    PRIVATE
    Catch2::Catch2WithMain
    moire_parallel_support
)

add_executable(transmission_incremental_tests tests/transmission_incremental_tests.cpp)
target_link_libraries(transmission_incremental_tests
    PRIVATE
    Catch2::Catch2WithMain
    moire_parallel_support
)
target_sources(transmission_incremental_tests PRIVATE
    ../src/prob_any_missing.cpp
    ../src/combination_indices_generator.cpp
)

add_executable(ecoi_marginal_tests tests/ecoi_marginal_tests.cpp)
target_link_libraries(ecoi_marginal_tests
    PRIVATE
    Catch2::Catch2WithMain
    moire_parallel_support
)

add_executable(marginal_proposal_tests tests/marginal_proposal_tests.cpp)
target_link_libraries(marginal_proposal_tests
    PRIVATE
    Catch2::Catch2WithMain
    moire_parallel_support
)

# Benchmark executables
add_executable(multivector_benchmarks1 benchmarks/multivector_benchmark1.cpp)
add_executable(multivector_benchmarks2 benchmarks/multivector_benchmark2.cpp)
target_link_libraries(multivector_benchmarks1 PRIVATE moire_parallel_support)
target_link_libraries(multivector_benchmarks2 PRIVATE moire_parallel_support)

add_executable(multivector_regression_benchmarks benchmarks/multivector_regression_benchmark.cpp)
target_link_libraries(multivector_regression_benchmarks PRIVATE moire_parallel_support)

add_executable(prob_any_missing_benchmarks benchmarks/prob_any_missing_benchmark.cpp)
target_link_libraries(prob_any_missing_benchmarks PRIVATE moire_parallel_support)
target_sources(prob_any_missing_benchmarks PRIVATE
    ../src/prob_any_missing.cpp
    ../src/combination_indices_generator.cpp
)

add_executable(pam_cache_calibration benchmarks/pam_cache_calibration.cpp)
target_link_libraries(pam_cache_calibration PRIVATE moire_parallel_support)
target_sources(pam_cache_calibration PRIVATE
    ../src/prob_any_missing.cpp
    ../src/combination_indices_generator.cpp
)

add_executable(ragged_multivector_benchmarks benchmarks/ragged_multivector_benchmark.cpp)
target_link_libraries(ragged_multivector_benchmarks PRIVATE moire_parallel_support)

add_executable(fused_mcmc_benchmarks benchmarks/fused_mcmc_benchmark.cpp)
target_link_libraries(fused_mcmc_benchmarks PRIVATE moire_parallel_support)

# Legacy test executables removed - using new Catch2-based framework

# Enable testing
enable_testing()

# Add tests
add_test(NAME multivector_tests COMMAND multivector_tests)
add_test(NAME distance_matrix_tests COMMAND distance_matrix_tests)
add_test(NAME observation_model_tests COMMAND observation_model_tests)
add_test(NAME prob_any_missing_tests COMMAND prob_any_missing_tests)
add_test(NAME zeta_mobius_tests COMMAND zeta_mobius_tests)
add_test(NAME transmission_incremental_tests COMMAND transmission_incremental_tests)
add_test(NAME ecoi_marginal_tests COMMAND ecoi_marginal_tests)
add_test(NAME marginal_proposal_tests COMMAND marginal_proposal_tests)

# Set test properties
set_tests_properties(multivector_tests PROPERTIES
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

set_tests_properties(distance_matrix_tests PROPERTIES
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

set_tests_properties(zeta_mobius_tests PROPERTIES
    TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Custom targets
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS multivector_tests distance_matrix_tests observation_model_tests prob_any_missing_tests zeta_mobius_tests transmission_incremental_tests ecoi_marginal_tests marginal_proposal_tests
    COMMENT "Running all tests"
)

add_custom_target(run_benchmarks
    COMMAND ./multivector_benchmarks1
    COMMAND ./multivector_benchmarks2 100 100 100
    COMMAND ./multivector_regression_benchmarks
    # Provide sensible defaults for argumented benchmarks
    COMMAND ./prob_any_missing_benchmarks 0.1 0.2 0.3 10
    COMMAND ./ragged_multivector_benchmarks 10 3 5,10,15
    DEPENDS multivector_benchmarks1 multivector_benchmarks2 multivector_regression_benchmarks prob_any_missing_benchmarks ragged_multivector_benchmarks
    COMMENT "Running all benchmarks (use MOIRE_BENCHMARK_CSV=1 for regression CSV)"
)

# Legacy test target removed - using new Catch2-based framework

# Coverage support
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(multivector_tests PRIVATE --coverage)
        target_compile_options(distance_matrix_tests PRIVATE --coverage)
        target_compile_options(zeta_mobius_tests PRIVATE --coverage)
        target_link_libraries(multivector_tests PRIVATE --coverage)
        target_link_libraries(distance_matrix_tests PRIVATE --coverage)
        target_link_libraries(zeta_mobius_tests PRIVATE --coverage)
    endif()
endif()
