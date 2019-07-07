cmake_minimum_required(VERSION 3.10.2)
project(OptimalTransport)

set(OptimalTransport_LIBRARIES OptimalTransport)

#Eigen 3 is required
if(NOT ITK_USE_SYSTEM_EIGEN)
  #Set Eigen3_DIR to the internal ITKEigen3 module.
  set(_internal_cmake_eigen3)
  list(GET ITKEigen3_INCLUDE_DIRS 0 _internal_cmake_eigen3)
  set(Eigen3_DIR "${_internal_cmake_eigen3}/itkeigen")
endif()
  find_package(Eigen3 REQUIRED)
  set(${PROJECT_NAME}_EXPORT_CODE_INSTALL
"${${PROJECT_NAME}_EXPORT_CODE_INSTALL} set(Eigen3_DIR \"${Eigen3_DIR}\") find_package(Eigen3 REQUIRED CONFIG)")


#Fetch optimal transport code
include(FetchContent)

set(OT_GIT_REPOSITORY "https://bitbucket.org/suppechasper/optimaltransport.git")
FetchContent_Declare(
    ot_fetch
    GIT_REPOSITORY ${OT_GIT_REPOSITORY} )

if(NOT ot_fetch_POPULATED)
  FetchContent_Populate(ot_fetch)
  add_subdirectory( ${ot_fetch_SOURCE_DIR} ${ot_fetch_BINARY_DIR} )
endif()

#Need manual icludes for now: Fix ot CMakelists
INCLUDE_DIRECTORIES( ${ot_fetch_SOURCE_DIR}/.
                     ${ot_fetch_SOURCE_DIR}/lib
                     ${ot_fetch_SOURCE_DIR}/external/gmra/lib
                     ${ot_fetch_SOURCE_DIR}/external/gmra/external/utils/lib
                     ${ot_fetch_SOURCE_DIR}/external/gmra/external/kmeans/lib
                     ${ot_fetch_SOURCE_DIR}/external/gmra/external/randomsvd/lib
                     ${ot_fetch_SOURCE_DIR}/external/gmra/external/metrics/lib
                     ${EIGEN3_INCLUDE_DIR} 
                     )


if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED)
  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
  include(ITKModuleExternal)
else()
  set(ITK_DIR ${CMAKE_BINARY_DIR})
  itk_module_impl()
endif()


