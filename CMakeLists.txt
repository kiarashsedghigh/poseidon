cmake_minimum_required(VERSION 3.20)
project(poseidon3_examples LANGUAGES C)

option(SD_PERMUTE "Compile SD with the [2,0,1] rotation" ON)
option(SD_SCALE "Compile SD with multiplier scaling" ON)

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

foreach(target IN ITEMS poseidon3-sd perf)
    target_compile_definitions(${target} PRIVATE
        SD_PERMUTE=$<BOOL:${SD_PERMUTE}> SD_SCALE=$<BOOL:${SD_SCALE}>)
endforeach()

# General-width experimental profiles with the existing StarkWare field.
add_executable(poseidon-profiles examples/poseidon_profiles.c sources/poseidon_profiles.c ${POSEIDON_FIELD_SOURCES})
target_include_directories(poseidon-profiles PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/sources")
target_compile_features(poseidon-profiles PRIVATE c_std_99)
target_compile_definitions(poseidon-profiles PRIVATE SD_PERMUTE=$<BOOL:${SD_PERMUTE}> SD_SCALE=$<BOOL:${SD_SCALE}>)
if(POSEIDON_ISO_C)
    target_compile_definitions(poseidon-profiles PRIVATE ISO_C=1)
endif()

add_executable(profiles-perf bench/profiles_perf.c sources/poseidon_profiles.c ${POSEIDON_FIELD_SOURCES})
target_include_directories(profiles-perf PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/sources")
target_compile_features(profiles-perf PRIVATE c_std_99)
target_compile_definitions(profiles-perf PRIVATE SD_PERMUTE=$<BOOL:${SD_PERMUTE}> SD_SCALE=$<BOOL:${SD_SCALE}>)
if(POSEIDON_ISO_C)
    target_compile_definitions(profiles-perf PRIVATE ISO_C=1)
endif()

add_executable(perf-80-t3 bench/perf_80_t3.c sources/poseidon_profiles.c ${POSEIDON_FIELD_SOURCES})
target_include_directories(perf-80-t3 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/sources")
target_compile_features(perf-80-t3 PRIVATE c_std_99)
target_compile_definitions(perf-80-t3 PRIVATE SD_PERMUTE=$<BOOL:${SD_PERMUTE}> SD_SCALE=$<BOOL:${SD_SCALE}>)
if(POSEIDON_ISO_C)
    target_compile_definitions(perf-80-t3 PRIVATE ISO_C=1)
endif()
