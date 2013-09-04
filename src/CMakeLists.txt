cmake_minimum_required (VERSION 2.6)
project (bmrfc)
set (bmrfcpp_VERSION_MAJOR 1)
set (bmrfcpp_VERSION_MINOR 0)

find_package(Boost COMPONENTS system thread REQUIRED)
find_package(HDF5 REQUIRED)
include(CMakeForceCompiler)
##find_package(OpenMP REQUIRED)


###### Compiler Optimizations
set(INTEL $ENV{INTEL_LICENSE_FILE})
if(NOT INTEL)
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -O3")
	# For Ba
	#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -O3 -march=opteron-sse3")
	# For profiling / debugging
	#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -O3 -g") # -pg 
else()
	MESSAGE(STATUS "** Using Intel Compiler **")
	CMAKE_FORCE_C_COMPILER(icc "Intel C Compiler")
	CMAKE_FORCE_CXX_COMPILER(icpc "Intel C++ Compiler")
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -O0")
endif()

####### Include Directories
include_directories(
  ${HDF5_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/libbmrfc
)

  
####### LibBMRFC
ADD_LIBRARY(bmrfc
  libbmrfc/BmrfHostAlgorithm.cpp
  libbmrfc/MatlabFile.cpp
  libbmrfc/IoFactory.cpp
  libbmrfc/MatlabFile.cpp
  libbmrfc/Timer.cpp
  libbmrfc/bmrf_functions.cpp
  libbmrfc/cdf_functions.cpp
  libbmrfc/matvec_functions.cpp
  libbmrfc/set_functions.cpp
  )
if(INTEL)
    set_target_properties(bmrfc PROPERTIES COMPILE_FLAGS "-openmp")
elseif(OPENMP_FOUND)
    set_target_properties(bmrfc PROPERTIES COMPILE_FLAGS "-fopenmp")
endif()
target_link_libraries(bmrfc
  log4cxx
  ${HDF5_LIBRARIES}
  ${Boost_THREAD_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  )

###### C Application Build
add_executable(bmrf-cpu bmrfcpu/bmrfcpu-main.cpp)
target_link_libraries(bmrf-cpu
  bmrfc
  )
if(INTEL)
	set_target_properties(bmrf-cpu PROPERTIES LINK_FLAGS " -openmp")
elseif(OPENMP_FOUND)
	set_target_properties(bmrf-cpu PROPERTIES LINK_FLAGS " -fopenmp")
endif()


####### CPU Tests
add_executable(bmrfc-tests
   bmrfc-tests/bmrfc-tests.cpp
   bmrfc-tests/cputests.cpp
   )
target_link_libraries(bmrfc-tests
   bmrfc
   )
if(INTEL)
	set_target_properties(bmrfc-tests PROPERTIES LINK_FLAGS " -openmp")
elseif(OPENMP_FOUND)
	set_target_properties(bmrfc-tests PROPERTIES LINK_FLAGS " -fopenmp")
endif()
