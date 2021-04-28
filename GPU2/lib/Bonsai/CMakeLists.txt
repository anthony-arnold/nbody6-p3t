CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

option(CUDA_VERBOSE_PTXAS
  "On to enable verbose output from the PTXAS assembler."
  OFF
  )

option(CUDA_DEVICE_DEBUGGING
  "On to enable debugging device code."
  OFF
  )

option(USE_MPI
  "On to build support for MPI"
  ON
  )

option(USE_MPIMT
  "On to build support for multi-threaded MPI"
  ON
)

option(USE_GALACTICS
  "Include John Dubinsky galatics IC generator"
  OFF
  )

option(USE_GALACTICS_IFORT
  "Set to ON if galactics is compipled with ifort"
  OFF
)

option(USE_THRUST
  "On to use thrust for sorting"
  OFF
  )

option(USE_CUB
    "On to use CUB sorting"
  ON
  )

option(COMPILE_SM35
  "On to enable SM35 compilation"
  ON
  )

option(USE_DUST
  "On to enable dust rendering"
  OFF
  )

option(CUDA_KEEP_INTERMEDIATE_FILES
  "On to enable -keep"
  OFF
  )

FIND_PACKAGE(CUDA REQUIRED)
if (USE_CUB)
  FIND_PACKAGE(CUB)
endif()

add_definitions(-std=c++11)
set(EXTRA_MPI_LINK_FLAGS)
if (USE_MPI)
  add_definitions(-DUSE_MPI)
  add_definitions(-fPIE -Wno-main)
  set(EXTRA_MPI_LINK_FLAGS   -pie -rdynamic)
if (USE_MPIMT)
  add_definitions(-D_MPIMT)
endif (USE_MPIMT)
endif (USE_MPI)

if (CUB_FOUND)
  add_definitions(-DUSE_CUB)
  include_directories(${CUB_INCLUDE_DIRS})
else ()
  add_definitions(-DUSE_THRUST)
endif ()
include_directories(${CUDA_TOOLKIT_ROOT_DIR}/../../thrust)


if (USE_DUST)
  add_definitions(-DUSE_DUST)
endif (USE_DUST)
set(LIBRARY_NAME bonsai)

set(BINARY_CLR_NAME bonsai_clrshm)
set (CLR_CCFILES
  src/bonsai_clrshm.cpp
  )

set(BINARY_IO_NAME  bonsai_io)
set (IO_CCFILES
  src/bonsai_io.cpp
  src/anyoption.cpp
  )

set(CMAKE_DEBUG_POSTFIX "D")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set (CCFILES
  src/build.cpp
  src/compute_properties.cpp
  src/gpu_iterate.cpp
  src/libraryInterface.cpp
  src/load_kernels.cpp
  src/main.cpp
  src/anyoption.cpp
  src/octree.cpp
  src/parallel.cpp
  src/sort_bodies_gpu.cpp
  src/log.cpp
  src/hostConstruction.cpp
  src/tipsyIO.cpp
  )

set (HFILES
  include/my_cuda_rt.h
  include/my_ocl.h
  include/node_specs.h
  include/octree.h
  include/build.h
  include/tipsydefs.h
  include/vector_math.h
  include/depthSort.h
  )

set (CUFILES
  CUDAkernels/build_tree.cu
  CUDAkernels/compute_propertiesD.cu
  CUDAkernels/parallel.cu
  CUDAkernels/sortKernels.cu
  CUDAkernels/scanKernels.cu
  CUDAkernels/timestep.cu
  CUDAkernels/depthSort.cu
  CUDAkernels/dev_direct_gravity.cu
  )

if (COMPILE_SM35)
  set (CUFILES
    ${CUFILES}
    CUDAkernels/dev_approximate_gravity_warp_new.cu
    )
  add_definitions("-DKEPLER")
  set(GENCODE
    -gencode arch=compute_35,code=sm_35
    -gencode arch=compute_35,code=compute_35
    -gencode arch=compute_50,code=compute_50
    -gencode arch=compute_50,code=sm_50
    -ftz=true -lineinfo)
else (COMPILE_SM35)
  set (CUFILES
    ${CUFILES}
	CUDAkernels/dev_approximate_gravity_warp_fermi.cu
	)
  set(GENCODE -gencode arch=compute_20,code=sm_20 -gencode arch=compute_20,code=compute_20)
endif (COMPILE_SM35)

set (CUHFILES
  CUDAkernels/support_kernels.cu
  CUDAkernels/scanKernels.cu
  )

set (PROFFILES
  # profiling/cuxTimer.cu
  # profiling/cuxTimer_host.cu
  profiling/cuxTimer.h
  profiling/cuxTimer_private.h
  profiling/warp_functions.h
  profiling/derived_atomic_functions.h
  profiling/bonsai_timing.h
  )

source_group("CUDA Source Files" FILES ${CUFILES} ${CUHFILES})
source_group("CUDA Include Files" FILES ${CUHFILES})
source_group("CUDA Profiling Tools" FILES ${PROFFILES})

if (CUDA_KEEP_INTERMEDIATE_FILES)
  set(KEEP -keep)
endif(CUDA_KEEP_INTERMEDIATE_FILES)

if (CUDA_VERBOSE_PTXAS)
  set(VERBOSE_PTXAS --ptxas-options=-v)
endif (CUDA_VERBOSE_PTXAS)

if (CUDA_DEVICE_DEBUGGING)
  set(DEVICE_DEBUGGING -G)
endif (CUDA_DEVICE_DEBUGGING)

if (WIN32)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
endif (WIN32)

add_definitions(-DHUGE=1.0e20)

#The SSE support for the local tree walks
if(WIN32)
  add_definitions( "/arch:SSE4")
 else (WIN32)

	#For PPC we need different vector flags
	if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
	    #Set compiler flags for IBM vector instructions
	elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
	    #This should support sse instructions
	    add_definitions( "-msse4  -march=native")
	else()
	    message(FATAL_ERROR "Unknown processor:" ${CMAKE_SYSTEM_PROCESSOR})
	endif()
endif (WIN32)

#if (USE_MPI)
	#The OpenMP Library and compiler flags
	FIND_PACKAGE(OpenMP REQUIRED)
	if(OPENMP_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
	#endif(USE_MPI)

set(lib_sources ${CCFILES})
list(REMOVE_ITEM lib_sources src/main.cpp)
cuda_add_library(${LIBRARY_NAME}
  ${lib_sources}
  ${HFILES}
  ${CUFILES}
  ${PROFFILES}
  OPTIONS ${GENCODE} ${VERBOSE_PTXAS} ${DEVICE_DEBUGGING} ${KEEP} -Xcompiler="-fPIC" -std=c++11
  )

if (USE_MPI)
    add_executable(${BINARY_CLR_NAME}
      ${CLR_CCFILES}
      )
    add_executable(${BINARY_IO_NAME}
      ${IO_CCFILES}
      )

    add_executable(bonsai_driver
      src/driver.cpp
      )
endif(USE_MPI)

if (USE_MPI)
    target_link_libraries(${BINARY_CLR_NAME}  -lrt ${EXTRA_MPI_LINK_FLAGS})
    target_link_libraries(${BINARY_IO_NAME}   -lrt ${EXTRA_MPI_LINK_FLAGS})
    target_link_libraries(bonsai_driver -ldl  -lrt ${EXTRA_MPI_LINK_FLAGS})
endif(USE_MPI)
