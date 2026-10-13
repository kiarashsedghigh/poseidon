cmake_minimum_required(VERSION 3.20)
project(poseidon3_examples LANGUAGES C)

option(POSEIDON_ISO_C "Use the portable ISO C field backend" OFF)

# A source list, not an additional build target.
set(POSEIDON_FIELD_SOURCES
    sources/f251.c
    sources/f251_int128.c
    sources/poseidon_rc.c
)

add_executable(poseidon3-target
    examples/poseidon3_pow5_standard.c
    sources/poseidon3_p5_plain.c
    ${POSEIDON_FIELD_SOURCES}
)

add_executable(poseidon3-sd
    examples/poseidon3_sd.c
    sources/posedion3_p5_sd.c
    ${POSEIDON_FIELD_SOURCES}
)

add_executable(perf
    bench/perf.c
    sources/poseidon3_p5_plain.c
    sources/posedion3_p5_sd.c
    ${POSEIDON_FIELD_SOURCES}
)

foreach(target IN ITEMS poseidon3-target poseidon3-sd perf)
    target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/sources")
    target_compile_features(${target} PRIVATE c_std_99)
    if(POSEIDON_ISO_C)
        target_compile_definitions(${target} PRIVATE ISO_C=1)
    endif()
endforeach()
