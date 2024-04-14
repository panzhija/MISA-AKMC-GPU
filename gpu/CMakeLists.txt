# set binary path and lib storage path for compiled binary.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_library(${KMC_GPU_LIB_NAME} STATIC
        gpu_simulate.hip
        gpu_simulate.h
        gpuError.hip
        gpuError.h
        DeviceVacRatesSolver.h
        )

target_compile_options(${KMC_GPU_LIB_NAME} PRIVATE -D__HIP_PLATFORM_HCC__= -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_HIP  -D__HIP__ -D__HIP_PLATFORM_AMD__= -I"${ROCM_PATH}/hip/include" -I"${ROCM_PATH}/llvm/bin/../lib/clang/13.0.0" -I${ROCM_PATH}/hsa/include)

target_include_directories(${KMC_GPU_LIB_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
        )

target_link_libraries(${KMC_GPU_LIB_NAME}
        ${PKMC_LIB_NAME}
        )

# install bin and libs
install(TARGETS ${KMC_GPU_LIB_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

# install header.
install(DIRECTORY "./"
        DESTINATION "include"
        FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
        )
