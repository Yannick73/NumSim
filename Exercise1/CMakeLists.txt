cmake_minimum_required(VERSION 3.8)

set(${PROJECT_NAME} "NumSim_Exc1")
# Define the project name.
set(${PROJECT_SOURCE_DIR} "src")
# Define the src directory

project(${PROJECT_NAME})

# Specify the name of the executable (${PROJECT_NAME} which is equal to what was set in the project() command).

# Also specify the source files.

add_executable(${PROJECT_NAME}

   main.cpp

)

# Add the project directory to include directories, to be able to include all project header files from anywhere

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})



# install numsim executable in build directory

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)

# Add additonial compile options to enable more warnings

add_compile_options(-Wall -Wextra)


# Set the version of the C++ standard to use, we use C++14, published in 2014

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_STANDARD_REQUIRED ON)


message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
