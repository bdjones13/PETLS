cmake_minimum_required(VERSION 3.16.3)

# change PETLS_INSTALL_PREFIX if you
# installed using "-DCMAKE_INSTALL_PREFIX=...", 
# e.g. /mnt/home/jones657/customlib
set(PETLS_INSTALL_PREFIX "/usr/local/") 


# include library and headers
set(PETLS_LIBRARIES ${PETLS_INSTALL_PREFIX}/lib/libpetls.a)
set(PETLS_INCLUDE_DIRS ${PETLS_INSTALL_PREFIX}/include/
                    ${PETLS_INSTALL_PREFIX}/include/petls_headers
                    ${PETLS_INSTALL_PREFIX}/include/petls_headers/core
                    ${PETLS_INSTALL_PREFIX}/include/petls_headers/eigenvalues
                    ${PETLS_INSTALL_PREFIX}/include/petls_headers/variants
                    )


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(CMAKE_CXX_FLAGS "-O3 -march=native -mtune=native -ffast-math -Wall -Wextra -pthread -frounding-math ${COMPILER_IGNORE_WARNINGS}") # -ffast-math  bad for bdcsvd?

project(main LANGUAGES C CXX)

# get Eigen - this is required to link against Eigen, which is part of the library interface
# Note: commonly installed versions of eigen are not sufficient. Eigen 3.4 or above is required
#       Eigen 3.3 is missing critical features. 
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

# make the executable
add_executable(main main.cpp)
target_include_directories(main PUBLIC ${PETLS_INCLUDE_DIRS})
target_link_libraries(main ${Eigen_LIBRARIES} Eigen3::Eigen  ${PETLS_LIBRARIES})



