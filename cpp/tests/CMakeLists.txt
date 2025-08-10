# on HPCC to install library locall use the cmake command
# cmake -DCMAKE_INSTALL_PREFIX=/mnt/home/jones657/customlib .. && make all install

# make sure the following modules are loaded:
# GCCcore/11.3.0, CMake/3.23.1, Eigen/3.4.0

# to dynamically link in the library, there is another environment variable you may need to export
# $LD_LIBRARY_PATH
cmake_minimum_required(VERSION 3.16.3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

project(test_cpp LANGUAGES C CXX)

add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

link_libraries(${Eigen_LIBRARIES} petls)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/" "${CMAKE_CURRENT_SOURCE_DIR}/src/core" "${CMAKE_CURRENT_SOURCE_DIR}/src/variants")

include(CTest)
enable_testing()

# define function for adding many tests
function (pl_add_test testname) # optional extra argument for subfolder of src
  if("${ARGV1}" STREQUAL "")
    set(filename "./src/test_${testname}.cpp")
  else()
    set(filename "./src/${ARGV1}/test_${testname}.cpp")
  endif()
  message("filename=${filename}")
  add_executable("test_${testname}" ${filename})
  # target_precompile_headers("test_${testname}" REUSE_FROM petls)

  add_test( NAME "cpp-${testname}" 
            COMMAND "test_${testname}"
            WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/")
endfunction()

# # https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html

# the tests 
# pl_add_test(cpp core)
pl_add_test(base core)
pl_add_test(alpha variants)
# pl_add_test(psl variants)
pl_add_test(rips variants)
pl_add_test(flag variants)
pl_add_test(flipped variants)
pl_add_test(reduction core)
# pl_add_test(eigenvectors core)
pl_add_test(eigs_algs core)
pl_add_test(psl_new variants)


