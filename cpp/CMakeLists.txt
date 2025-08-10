cmake_minimum_required(VERSION 3.16.3)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 


option(PETLS_USE_LAPACK "Use LAPACK/OpenBLAS for linear algebra" ON)

if(PETLS_USE_LAPACK)
  find_package(OpenBLAS)
  if(OpenBLAS_FOUND)
    message(STATUS "found openblas and lapack")
    SET(LAPACK_FLAGS "-lopenblas -llapacke")
    add_definitions (-DEIGEN_USE_LAPACKE)
  endif()
endif()

set(COMPILER_IGNORE_WARNINGS "-Wno-unused-variable -Wno-unused-parameter -Wno-unused-private-field -Wno-maybe-uninitialized") # most of these are necessary to repress warnings from Eigen and other dependencies

set(CMAKE_CXX_FLAGS "${LAPACK_FLAGS} -fPIC -O3 -march=native -mtune=native -ffast-math -Wall -Wextra -pthread -frounding-math ${COMPILER_IGNORE_WARNINGS}") # removed -Werror for Gudhi

# pthread from flagser
# set(CMAKE_BUILD_TYPE DEBUG) # CGAL warns to use release for performance reasons
add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

project(petls VERSION 0.1.0 LANGUAGES C CXX)


set(Boost_NO_BOOST_CMAKE ON) # on some manylinux systems, BoostConfig.cmake is not found, this is a workaround
find_package(Boost REQUIRED COMPONENTS program_options thread)


# Include dir
include_directories(/usr/local/include)

# Src
AUX_SOURCE_DIRECTORY(src SRC_FILES)

# get eigen
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

if (PETLS_USE_ALPHA_COMPLEX)
  include(alpha.cmake)
endif()

include(gudhi.cmake)


#### Begin Directed Flag Complex ####
find_package(Threads REQUIRED)
##### End Directed Flag Complex #####


# Headers
set(PROJECT_SOURCE_DIR  ${CMAKE_CURRENT_LIST_DIR}/src)
set(PROJECT_INCLUDE_DIR  ${CMAKE_CURRENT_LIST_DIR}/include/petls_headers)
set(petls_INCLUDE_DIR ${PROJECT_INCLUDE_DIR})
# Source files

set(SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/petls.cpp
    ${PROJECT_SOURCE_DIR}/core/Complex.cpp
    ${PROJECT_SOURCE_DIR}/eigenvalues/eigs_algorithms.cpp

    ${PROJECT_SOURCE_DIR}/variants/Alpha.cpp
)

set(HEADER_FILES 
    ${PROJECT_INCLUDE_DIR}/petls.hpp
    ${PROJECT_INCLUDE_DIR}/core/FilteredBoundaryMatrix.hpp
    ${PROJECT_INCLUDE_DIR}/core/Complex.hpp
    ${PROJECT_INCLUDE_DIR}/core/up_algorithms.hpp
    
    ${PROJECT_INCLUDE_DIR}/eigenvalues/eigs_algorithms.hpp

    ${PROJECT_INCLUDE_DIR}/variants/Alpha.hpp # Alpha Complex depends on CGAL
    ${PROJECT_INCLUDE_DIR}/variants/dFlag.hpp
    ${PROJECT_INCLUDE_DIR}/variants/Rips.hpp
    
    ${PROJECT_INCLUDE_DIR}/variants/newPSL.hpp
    ${PROJECT_INCLUDE_DIR}/variants/PersistentSheafLaplacian.hpp
    ${PROJECT_INCLUDE_DIR}/variants/sheaf_simplex_tree.hpp
)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES})

# Add library
add_library(petls "")
target_sources(petls
  PRIVATE ${SOURCE_FILES}
  PUBLIC ${HEADER_FILES}

)
set_target_properties(petls PROPERTIES LINKER_LANGUAGE CXX)


# Include directories
target_include_directories(petls PUBLIC  ${PROJECT_INCLUDE_DIR}
                                                        ${PROJECT_INCLUDE_DIR}/core
                                                        ${PROJECT_INCLUDE_DIR}/eigenvalues
                                                        ${PROJECT_INCLUDE_DIR}/variants)
##

target_link_libraries (petls ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})


if (TBB_FOUND)
  target_link_libraries(petls ${TBB_LIBRARIES})
endif(TBB_FOUND)
if (GMP_FOUND)
  target_link_libraries(petls ${GMP_LIBRARIES})
endif(GMP_FOUND)
if (GMPXX_FOUND)
  target_link_libraries(petls ${GMPXX_LIBRARIES})
endif(GMPXX_FOUND)
                                                        
target_link_libraries(petls Eigen3::Eigen)
 
# Install
install(TARGETS petls DESTINATION lib)

# # Install the headers
install(FILES include/petls DESTINATION include)

install(DIRECTORY include/petls_headers DESTINATION include)
install(DIRECTORY include DESTINATION include)

# can force debug with cmake -S .. -B . -D CMAKE_BUILD_TYPE=Debug
if(CMAKE_BUILD_TYPE MATCHES Debug)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)     
endif()

