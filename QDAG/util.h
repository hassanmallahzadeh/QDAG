//
// Created by stefan on 04.12.19.
//

#ifndef DD_PACKAGE_UTIL_H
#define DD_PACKAGE_UTIL_H
#include "IIC-JKU/DDpackage.h"
#include <string>
#include <cassert>

// X gate matrix
constexpr dd::Matrix2x2 Xmat = {{{ 0, 0 }, { 1, 0 } }, {{ 1, 0 }, { 0, 0 } }};
// Z gate matrix
constexpr dd::Matrix2x2 Zmat = {{{ 1, 0 }, { 0, 0 } }, {{ 0, 0 }, { -1, 0 } }};
// Y gate matrix
constexpr dd::Matrix2x2 Ymat = {{{ 0, 0 }, { 0, -1 } }, {{ 0, 1 }, { 0, 0 } }};
// Hadamard gate matrix
constexpr dd::Matrix2x2 Hmat = {{{ dd::SQRT_2, 0 }, { dd::SQRT_2,  0 }},
                            {{ dd::SQRT_2, 0 }, { -dd::SQRT_2, 0 }}};
// all 1 gate matrix (used for test)
constexpr dd::Matrix2x2 oneMat = {{{ 1, 0 }, { 1, 0 }},
                                {{ 1, 0 }, { 1, 0 }}};
#endif //DD_PACKAGE_UTIL_H
