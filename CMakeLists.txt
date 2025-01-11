cmake_minimum_required(VERSION 3.10)

project(MeshRenderer)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(VTK REQUIRED COMPONENTS
  vtkCommonCore
  vtkCommonDataModel
  vtkIOGeometry
  vtkRenderingCore
  vtkInteractionStyle
  vtkRenderingOpenGL2
)

include(${VTK_USE_FILE})

add_executable(MeshRenderer main.cpp)

target_link_libraries(MeshRenderer PRIVATE ${VTK_LIBRARIES})