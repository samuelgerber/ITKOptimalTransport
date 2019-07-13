cmake_minimum_required(VERSION 3.10.2)
project(OptimalTransport)

#set(OptimalTransport_LIBRARIES OptimalTransport)



#Eigen 3 is required
if(NOT ITK_USE_SYSTEM_EIGEN)
  #Set Eigen3_DIR to the internal ITKEigen3 module.
  set(_internal_cmake_eigen3)
  list(GET ITKEigen3_INCLUDE_DIRS 0 _internal_cmake_eigen3)
  set(Eigen3_DIR "${_internal_cmake_eigen3}/itkeigen")
endif()
  find_package(Eigen3 REQUIRED)
  #find_package (Eigen3 3.3 REQUIRED NO_MODULE)
  set(${PROJECT_NAME}_EXPORT_CODE_INSTALL
"${${PROJECT_NAME}_EXPORT_CODE_INSTALL} set(Eigen3_DIR \"${Eigen3_DIR}\") find_package(Eigen3 REQUIRED CONFIG)")



if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED)
  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
  include(ITKModuleExternal)
else()
  set(ITK_DIR ${CMAKE_BINARY_DIR})
  itk_module_impl()
endif()



add_subdirectory("examples")
