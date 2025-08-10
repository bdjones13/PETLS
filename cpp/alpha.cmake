# put all of the logic for including Alpha Laplacian
# need to check for CGAL, other dependencies, install Gudhi


#### GUDHI Alpha complex requires CGAL >= 4.7.0 ####
find_package(CGAL QUIET)

# Requires CGAL versions > 4.11
if (CGAL_FOUND)



  #### GUDHI requires Boost - Alpha complex requires Boost program_options thread ####
  # find_package(Boost 1.78.0 REQUIRED COMPONENTS program_options thread HINTS $ENV{BOOST_ROOT})

  # BOOST ISSUE result_of vs C++11
  add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)
  # BOOST ISSUE with Libraries name resolution under Windows
  add_definitions(-DBOOST_ALL_NO_LIB)
  # problem with Visual Studio link on Boost program_options
  add_definitions( -DBOOST_ALL_DYN_LINK )
  # problem on Mac with boost_system and boost_thread
  add_definitions( -DBOOST_SYSTEM_NO_DEPRECATED )
  include_directories(${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS})


  if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
    # add sources for Alpha  
    include( ${CGAL_USE_FILE} )
    #### Optional GMP and GMPXX for CGAL ####
    find_package(GMP REQUIRED)
    if(GMP_FOUND)
      include_directories(${GMP_INCLUDE_DIR})
      find_package(GMPXX REQUIRED)
    endif()
    if(GMPXX_FOUND)
      include_directories(${GMPXX_INCLUDE_DIR})
    endif()
    message("++ Found CGAL, GMP, GMPXX - will include the Alpha.")
    # set(ALPHA_SRC "src/variants/Alpha.cpp")
    # set(ALPHA_INCLUDE "include/petls_headers/variants/Alpha.hpp")
    # set(ALPHA_PYTHON "src/variants/Alpha.cpp")
    add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK) # Allows Valgrind to run with CGAL/Gudhi/Alpha complex
 
    add_definitions(-DPETLS_USE_ALPHA_COMPLEX)
    set_property(GLOBAL PROPERTY PETLS_USING_ALPHA ON)
    #### Optional TBB for CGAL and GUDHI ####
    set(TBB_FIND_QUIETLY ON)
    find_package(TBB)
    if (TBB_FOUND)
      include(${TBB_USE_FILE})
    endif()

  else()
    message("++ CGAL Version is ${CGAL_VERSION}, but the Alpha requires version > 4.11.")
    message("++   The Alpha class will not be installed.")
    # set(ALPHA_PYTHON "src/variants/Alpha_placeholder.cpp") # Caution: Dummy class wraps Simplex
  endif()
else()
  message("++ CGAL version > 4.11 is required for the Alpha. The Alpha class will not be installed.")
  # set(ALPHA_PYTHON "src/variants/Alpha_placeholder.cpp") # Caution: Dummy class wraps Simplex
endif()

