cmake_minimum_required( VERSION 2.8 )

project( UtilLib )

# Find the mpi package
find_package(MPI REQUIRED)

# Use MPI compiler and linker
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

# Use mpicxx to replace 
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# Include header file directory to the search path
include_directories ( "${PROJECT_SOURCE_DIR}/include" )

add_library( UtilLib STATIC ${PROJECT_SOURCE_DIR}/src/SystemUtils.cpp
							${PROJECT_SOURCE_DIR}/src/IOUtils.cpp
							${PROJECT_SOURCE_DIR}/src/InterpUtils.cpp )

# Installation Instructions
# lib is interpreted as relative to CMAKE_INSTALL_PREFIX
FILE(GLOB headers "${PROJECT_SOURCE_DIR}/include/*.h")
install( FILES ${headers} DESTINATION include )
install( TARGETS UtilLib DESTINATION lib )
     
