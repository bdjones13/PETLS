# install gudhi
set(WITH_GUDHI_THIRD_PARTY OFF CACHE BOOL "docstring")
set(HERA_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/Gudhi_modified/ext/hera/include/" CACHE STRING "docstring" FORCE)
include(FetchContent)
FetchContent_Declare(Gudhi
                     SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/Gudhi_modified/" )


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include/Gudhi_modified/src/cmake/modules/")

FetchContent_GetProperties(Gudhi)

FetchContent_MakeAvailable(Gudhi)
include_directories("${gudhi_SOURCE_DIR}/include/" ${gudhi_BINARY_DIR} )
if(NOT gudhi_POPULATED)
    FetchContent_Populate(Gudhi)
    add_subdirectory(${gudhi_SOURCE_DIR} ${gudhi_BINARY_DIR})
    include_directories(${gudhi_SOURCE_DIR} ${gudhi_BINARY_DIR} )
endif()