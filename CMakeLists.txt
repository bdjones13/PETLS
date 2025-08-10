cmake_minimum_required(VERSION 3.16.3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
set(COMPILER_IGNORE_WARNINGS "-Wno-unused-variable -Wno-unused-parameter -Wno-unused-private-field")


option(PETLS_USE_ALPHA_COMPLEX "Enable Alpha complex functionality (requires Gudhi + CGAL)" ON)
# To set this to on, use "cmake -DPETLS_USE_ALPHA_COMPLEX=ON .."

if (PETLS_USE_ALPHA_COMPLEX MATCHES OFF)
    message(STATUS "PETLS_USE_ALPHA_COMPLEX IS SET TO OFF - Alpha complex will not be available in Python or C++")
endif()

project(petls)

# https://github.com/pypa/cibuildwheel/issues/639
#   Development.Embed not supported on manylinux so cannot only require "Development"
#   On local machine, requiring Development.Module only fails,
#   while just "Development" (which implies Development.Module and Development.Embed) works 
if(CMAKE_BUILD_TYPE MATCHES Debug)
    SET(Python_FIND_VIRTUALENV ONLY)
    find_package(Python COMPONENTS Interpreter Development REQUIRED ) 

else()
    SET(Python_FIND_VIRTUALENV FIRST)
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED )  
endif()

message("python found: ${Python_FOUND}, interpreter: ${Python_INTERPRETER}, library dirs: ${Python_LIBRARY_DIRS}, python_libraries: ${Python_LIBRARIES}")

add_subdirectory(cpp)

get_cmake_property(PETLS_USING_ALPHA PETLS_USING_ALPHA)
if (PETLS_USING_ALPHA) # set in cpp/alpha.cmake, will only be set if also PETLS_USE_ALPHA_COMPLEX is ON
    add_definitions(-DPETLS_USE_ALPHA_COMPLEX)
endif()
include(cpp/gudhi.cmake) # this is a fix to make sure gudhi_SOURCE_DIR and gudhi_BINARY_DIR are defined correctly at this scope
# endif()

include_directories("${CMAKE_SOURCE_DIR}/cpp/include/petls_headers")
include_directories("${CMAKE_SOURCE_DIR}/python")


set(PROJECT_SOURCE_DIR "cpp/src")
set(PROJECT_INCLUDE_DIR "cpp/include/petls_headers")
set(PYTHON_INCLUDE_DIR "src/")

set(SOURCE_FILES 
    ${PROJECT_SOURCE_DIR}/petls.cpp
    ${PROJECT_SOURCE_DIR}/core/Complex.cpp
    ${PROJECT_SOURCE_DIR}/eigenvalues/eigs_algorithms.cpp
    ${PROJECT_SOURCE_DIR}/variants/PersistentSheafLaplacian.cpp
    ${PROJECT_SOURCE_DIR}/variants/Alpha.cpp
)

set(HEADER_FILES 
    ${PROJECT_INCLUDE_DIR}/petls.hpp

    ${PROJECT_INCLUDE_DIR}/core/FilteredBoundaryMatrix.hpp
    ${PROJECT_INCLUDE_DIR}/core/Complex.hpp

    ${PROJECT_INCLUDE_DIR}/eigenvalues/eigs_algorithms.hpp
    ${PROJECT_INCLUDE_DIR}/core/up_algorithms.hpp
    ${PROJECT_INCLUDE_DIR}/variants/Alpha.hpp # Note - Alpha included above if CGAL + Gudhi installed
    ${PROJECT_INCLUDE_DIR}/variants/dFlag.hpp
    ${PROJECT_INCLUDE_DIR}/variants/Rips.hpp
    ${PROJECT_INCLUDE_DIR}/variants/newPSL.hpp
)


set(PYTHON_FILES 
    ${PYTHON_INCLUDE_DIR}/core/Complex.cpp
    ${PYTHON_INCLUDE_DIR}/variants/Alpha.cpp
    ${PYTHON_INCLUDE_DIR}/variants/dFlag.cpp
    ${PYTHON_INCLUDE_DIR}/variants/Rips.cpp
    ${PYTHON_INCLUDE_DIR}/variants/NewPersistentSheafLaplacian.cpp
    ${PYTHON_INCLUDE_DIR}/_petls.cpp    
)

# Set up such that XCode organizes the files
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )


include(pybind11.cmake) # supposedly find pybind11 after find python is best https://stackoverflow.com/questions/73053479/python-version-mismatch-even-though-cmake-reports-having-found-the-correct-versi

python_add_library(_petls MODULE
 	${SOURCE_FILES}
    ${HEADER_FILES}
    ${PYTHON_FILES}
)


target_link_libraries(_petls PUBLIC pybind11::headers Eigen3::Eigen ${GMP_LIBRARIES} ${GMPXX_LIBRARIES} petls)

target_include_directories(_petls PUBLIC "${persistentlaplacians_SOURCE_DIR}include")
target_include_directories(_petls PUBLIC "${gudhi_SOURCE_DIR}/include/" ${gudhi_BINARY_DIR} )

message("PYTHON_LIBRARY_DIRS: ${PYTHON_LIBRARY_DIRS}")
message("Python_EXECUTABLE:${Python_EXECUTABLE}")
if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(Python_EXECUTABLE "${CMAKE_SOURCE_DIR}/.venv/bin/python")
    set(PYTHON_LIBRARY_DIRS "${CMAKE_SOURCE_DIR}/.venv/lib/python3.12/site-packages")
    install(TARGETS _petls DESTINATION "${PYTHON_LIBRARY_DIRS}/petls") 
else() 
    install(TARGETS _petls DESTINATION petls)
endif()


include(CTest)
enable_testing()

# define function for adding many tests
function (pl_add_test_python testname) # optional extra argument for subfolder of src
  if("${ARGV1}" STREQUAL "")
    set(folder "${CMAKE_CURRENT_SOURCE_DIR}/tests/core/")
  else()
    set(folder "${CMAKE_CURRENT_SOURCE_DIR}/tests/${ARGV1}/")
  endif()
  
  set(filename "${folder}test_${testname}.py")
  message("filename=${filename}")
  add_test(NAME "python-${testname}" 
            COMMAND ${Python_EXECUTABLE} -m pytest ${filename}
            WORKING_DIRECTORY "${folder}")
endfunction()

pl_add_test_python("base")
pl_add_test_python("profile")
pl_add_test_python("rips" "variants")
pl_add_test_python("sst" "variants")
pl_add_test_python("psl" "variants")
