cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_FLAGS -std=c++14)
project(QDAG)
add_executable(
QDAG
QFT.cpp
IIC-JKU/DDcomplex.cpp
IIC-JKU/DDexport.cpp
IIC-JKU/DDpackage.cpp
QFT-DDgenerator.cpp
main.cpp
)