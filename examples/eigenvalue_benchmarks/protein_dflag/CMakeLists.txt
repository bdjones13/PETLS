# If you do not have root permissions, you may have to the install library locally use the cmake command
# cmake -DCMAKE_INSTALL_PREFIX=/mnt/home/username/customlib .. && make all install

# make sure the following modules are loaded:
# GCCcore/11.3.0, CMake/3.23.1, Eigen/3.4.0

cmake_minimum_required(VERSION 3.16.3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(COMPILER_IGNORE_WARNINGS "-Wno-unused-variable -Wno-unused-parameter -Wno-unused-private-field")
option(PETLS_USE_LAPACK "helptext" ON)

if(PETLS_USE_LAPACK)
  find_package(OpenBLAS)
  if(OpenBLAS_FOUND)
    message(STATUS "found openblas and lapack")
    SET(LAPACK_FLAGS "-lopenblas -llapacke")
    add_definitions (-DEIGEN_USE_LAPACKE)
  endif()
endif()
set(CMAKE_CXX_FLAGS "-g -O3 ${LAPACK_FLAGS} -march=native -mtune=native -ffast-math -Wall -Wextra -pthread -frounding-math ${COMPILER_IGNORE_WARNINGS}") # -ffast-math  bad for bdcsvd?


include(GNUInstallDirs)
include(FetchContent)
FetchContent_Declare(eigen 
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE )
SET(BUILD_TESTING OFF)
SET(EIGEN_BUILD_TESTING OFF)  
FetchContent_GetProperties(eigen)
FetchContent_MakeAvailable(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
    add_subdirectory(${Eigen_SOURCE_DIR} ${Eigen_BINARY_DIR})
    include_directories(${Eigen_SOURCE_DIR} ${Eigen_BINARY_DIR} )
endif()


project(cpp_eigs LANGUAGES C CXX)
set(PETLS_LINK /mnt/home/jones657/customlib/lib/libpetls.a)
link_libraries(${Eigen_LIBRARIES} ${PETLS_LINK})

include_directories("/mnt/home/jones657/customlib/include/petls_headers/core" "/mnt/home/jones657/customlib/include/petls_headers/variants" ${EIGEN_INCLUDE_DIRS})# "/mnt/home/jones657/customlib/include/petls_headers" "/mnt/home/jones657/customlib/lib")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../../cpp/include/spectra/include/")

add_executable(cpp_eigs cpp_eigs.cpp)
target_link_libraries(cpp_eigs Eigen3::Eigen)



