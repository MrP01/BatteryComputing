cmake_minimum_required(VERSION 3.16)
project(heatfun)

add_definitions("-std=c++20")
add_definitions("-O3")  # optimisation level 3
add_definitions("-g")  # gcc adds debug information (confer `man gcc`)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(Qt6 REQUIRED COMPONENTS Core)
find_package(Qt6 REQUIRED COMPONENTS Widgets)
find_package(Qt6 REQUIRED COMPONENTS Charts)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

# list(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}")
# find_package(pybind11 REQUIRED)
# pybind11_add_module(heatfun solver/module.cpp solver/Solver.cpp solver/TschebFun.cpp)
# target_link_libraries(heatfun PUBLIC muparserx)
# target_compile_definitions(heatfun PRIVATE SILENT=1)

add_executable(main
  SpectralSolver/main.cpp
  SpectralSolver/Interface.cpp
  SpectralSolver/Solver.cpp
  ../HeatFun/solver/Interface.cpp
  ../HeatFun/solver/Solver.cpp
  ../HeatFun/solver/TschebFun.cpp
)
target_link_libraries(main PUBLIC Qt6::Core Qt6::Widgets Qt6::Charts muparserx)
