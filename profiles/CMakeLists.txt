# set binary path and lib storage path for compiled binary.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/profile_config.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/profile_config.h"
)

set(P_KMC_EXECUTE_BIN_NAME pmisa-kmc)
set(P_KMP_PROFILE_LIB_NAME profile)

add_library(${P_KMP_PROFILE_LIB_NAME} STATIC
        config/config_parsing.cpp
        config/config_parsing.h
        config/config_values.cpp
        config/config_values.h
        config/lattice_types_string.cpp
        config/lattice_types_string.h
        io/lattice_dump.cpp
        io/lattice_dump.h
        device.h
        m_event_listener.h
        m_event_listener.cpp
        m_event_hook.cpp
        m_event_hook.h
        pkmc.h
        pkmc.cpp
        seed_relocate.h
        )


target_compile_options(${P_KMP_PROFILE_LIB_NAME} PRIVATE -D__HIP_PLATFORM_HCC__= -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_HIP  -D__HIP__ -D__HIP_PLATFORM_AMD__= -I"${ROCM_PATH}/hip/include" -I"${ROCM_PATH}/llvm/bin/../lib/clang/13.0.0" -I${ROCM_PATH}/hsa/include)

target_include_directories(${P_KMP_PROFILE_LIB_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
        )

target_link_libraries(${P_KMP_PROFILE_LIB_NAME} PUBLIC
        kiwi::kiwi
        yaml-cpp
        ${PKMC_LIB_NAME}
        ${KMC_EXTRA_LIBS}
        )

add_executable(${P_KMC_EXECUTE_BIN_NAME} profile.cpp)

target_compile_options(${P_KMC_EXECUTE_BIN_NAME} PRIVATE -D__HIP_PLATFORM_HCC__= -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_HIP  -D__HIP__ -D__HIP_PLATFORM_AMD__= -I"${ROCM_PATH}/hip/include" -I"${ROCM_PATH}/llvm/bin/../lib/clang/13.0.0" -I${ROCM_PATH}/hsa/include)

target_link_libraries(${P_KMC_EXECUTE_BIN_NAME} PRIVATE
        ${P_KMP_PROFILE_LIB_NAME}
        )

# install bin and libs
install(TARGETS ${P_KMC_EXECUTE_BIN_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )
