cmake_minimum_required(VERSION 3.20)
project(poseidon3_examples LANGUAGES C)

# The C backend automatically uses __int128 when supported by the compiler.
option(POSEIDON_ISO_C "Use the portable ISO C field backend" OFF)

add_library(poseidon3_field STATIC
    sources/f251.c
    sources/f251_int128.c
    sources/poseidon_rc.c
)
target_include_directories(poseidon3_field PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/sources")
target_compile_features(poseidon3_field PUBLIC c_std_99)
if(POSEIDON_ISO_C)
    target_compile_definitions(poseidon3_field PRIVATE ISO_C=1)
endif()

add_library(poseidon3_standard_impl STATIC sources/poseidon_pow5_standard.c)
target_link_libraries(poseidon3_standard_impl PUBLIC poseidon3_field)

add_library(poseidon3_sd_impl STATIC sources/poseidon_sd.c)
target_link_libraries(poseidon3_sd_impl PUBLIC poseidon3_field)

# The baseline example also checks its output against identity/restored SD.
add_executable(poseidon3-target examples/poseidon3_pow5_standard.c)
target_link_libraries(poseidon3-target PRIVATE poseidon3_standard_impl poseidon3_sd_impl)

add_executable(poseidon3-sd examples/poseidon3_sd.c)
target_link_libraries(poseidon3-sd PRIVATE poseidon3_sd_impl)

include(CTest)
if(BUILD_TESTING)
    add_test(NAME poseidon3-baseline-comparison COMMAND poseidon3-target)
    set_tests_properties(poseidon3-baseline-comparison PROPERTIES
        PASS_REGULAR_EXPRESSION "PASS: baseline matches identity SD and restored SD")
    add_test(NAME poseidon3-sd-default COMMAND poseidon3-sd)
    add_test(NAME poseidon3-sd-no-permutation COMMAND poseidon3-sd 0 3 5 7)
    set_tests_properties(poseidon3-sd-default poseidon3-sd-no-permutation PROPERTIES
        PASS_REGULAR_EXPRESSION "Python-style hash output is state\\[1\\]")
endif()
