cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_FLAGS -std=c++14)
project(QDAG)
add_executable(
QDAG
QFT.cpp
DDcomplex.cpp
DDexport.cpp
DDpackage.cpp
originalmain.cpp
QFT-DDgenerator.cpp
)