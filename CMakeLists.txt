cmake_minimum_required(VERSION 3.18)

project(quantlop LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(
    _core
    cpp/src/bindings.cpp
    cpp/src/pauliword.cpp
    cpp/src/matvec.cpp
    cpp/src/hamiltonian.cpp
)
target_include_directories(_core PRIVATE cpp/include)

install(TARGETS _core LIBRARY DESTINATION quantlop)
