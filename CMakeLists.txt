cmake_minimum_required(VERSION 3.30)
project(Physics_Sandbox CUDA CXX)

set(CMAKE_CUDA_STANDARD 20)
set(wxBUILD_SHARED OFF)

add_subdirectory(wx-3.2.6)
add_executable(Physics_Sandbox main.cu movement.cu bitmap.cpp shape.cpp utils.cu draw.cu)
add_executable(Physics_Sandbox_gui display_window.cpp movement.cu bitmap.cpp shape.cpp utils.cu draw.cu)

find_package(CUDAToolkit)

target_link_libraries(Physics_Sandbox wx::net wx::core wx::base)
target_link_libraries(Physics_Sandbox_gui wx::net wx::core wx::base)
target_link_libraries(Physics_Sandbox cublas nppc nppial)
target_link_libraries(Physics_Sandbox_gui cublas nppc nppial)

set_target_properties(Physics_Sandbox PROPERTIES
                      CUDA_SEPARABLE_COMPILATION ON)

