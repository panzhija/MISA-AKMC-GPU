cmake_minimum_required(VERSION 3.0)

project(misa-kmc)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_HIP_STANDARD 11)

include(version.cmake)
include(config.cmake)
include(dependency.cmake)

if (KMC_HIP_ENABLE_FLAG)
    enable_language(HIP)
    add_subdirectory(gpu)
endif ()

add_subdirectory(models)
add_subdirectory(src)
add_subdirectory(profiles)
add_subdirectory(profile_serial)
# tools
if (KMC_TOOLS_BUILD_FLAG)
    add_subdirectory(tools)
endif ()

#if (KMC_TOOLS_BUILD_ENABLE_FLAG)
#    add_subdirectory(tools)
#endif ()
################################
## test configure
################################
if (KMC_TEST_BUILD_ENABLE_FLAG)
    MESSAGE(STATUS "TEST is enabled")
    enable_testing()
    add_subdirectory(tests)
endif ()
