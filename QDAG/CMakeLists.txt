cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_FLAGS -std=c++17)
project(QDAG)
add_executable(
QDAG
QFT.cpp
IIC-JKU/DDcomplex.cpp
IIC-JKU/DDexport.cpp
IIC-JKU/DDpackage.cpp
QFT-DDgenerator.cpp
QFT-Measurement.cpp
RegisterFactory.cpp
shorutil.cpp
ProbabilisticPeriodFinder.cpp
OldProbabilisticPeriodFinder.cpp
SB-ProbabilisticPeriodFinder.cpp
Factorizer.cpp
main.cpp
)