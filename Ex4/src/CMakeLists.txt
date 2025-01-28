cmake_minimum_required(VERSION 3.8)

# Define the project name.
project(numsim_3d)

# Specify the name of the executable (${PROJECT_NAME} which is equal to what was set in the project() command).
# Also specify the source files.
add_executable(${PROJECT_NAME}

   settings.cpp
   storage/array2D.cpp
   storage/array3D.cpp
   storage/field_variable.cpp
   discretization/staggered_grid.cpp
   discretization/discretization.cpp
   discretization/donor_cell.cpp
   discretization/central_differences.cpp
   discretization/partition_information.cpp
   #discretization/partition.cpp
   discretization/async_partition.cpp
   boundary/boundary.cpp
   boundary/dirichlet.cpp
   #boundary/neighbour_boundary.cpp
   #boundary/async_neighbour_boundary.cpp
   #output_writer/output_writer.cpp
   #output_writer/output_writer_text.cpp
   #output_writer/output_writer_paraview.cpp
   #output_writer/output_writer_paraview_parallel.cpp
   pressure_solver/pressure_solver.cpp
   pressure_solver/gauss_seidel.cpp
   #pressure_solver/sor.cpp
   #pressure_solver/checkerboard.cpp
   #pressure_solver/cg.cpp
   computation.cpp
   main.cpp

   #testmain.cpp
   
)

# Add the project directory to include directories, to be able to include all project header files from anywhere
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

# Search for the external package "VTK"
find_package(VTK REQUIRED)

# Output, if VTK was found
message("If VTK was found on the system: VTK_FOUND: ${VTK_FOUND}")

# If an installation of vtk was found on the system
if (VTK_FOUND)
  include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
  target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
endif(VTK_FOUND)

# Search for MPI
find_package(MPI REQUIRED)

# Output, if MPI was found
message("If MPI was found on the system: MPI_FOUND: ${MPI_FOUND}")

include_directories(${MPI_INCLUDE_PATH})

target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)

  set_target_properties(${PROJECT_NAME} PROPERTIES

    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")

endif()


if(MPI_LINK_FLAGS)

  set_target_properties(${PROJECT_NAME} PROPERTIES

    LINK_FLAGS "${MPI_LINK_FLAGS}")

endif()

option(PROFILE "Add flags to profile the program with gprof." OFF)

if(PROFILE)

  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pg -g")

  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pg -g")

endif()

if(DT_STATISTICS)
  message("Set dt-statics mode")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDT_STATISTICS")
endif()

if(SOLVER_STATISTICS)
  message("Set pressure-solver-statics mode")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSOLVER_STATISTICS")
endif()

if(TIMER)
  message("Set timing mode")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIMER")
endif()

if(GEOMETRY)
  message("Set geometry out")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGEOMETRY")
endif()

# Wall and Wextra for warings & some for optimization
add_compile_options(-Wall -Wextra -O3 -Ofast -lto -march=native)

# install numsim executable in build directory
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)


# Set the version of the C++ standard to use, we use C++14, published in 2014
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
