//
// Created by stefan on 04.12.19.
//

#ifndef DD_PACKAGE_UTIL_H
#define DD_PACKAGE_UTIL_H
#include "IIC-JKU/DDpackage.h"

// X gate matrix
constexpr dd::Matrix2x2 Xmat = {{{ 0, 0 }, { 1, 0 } }, {{ 1, 0 }, { 0, 0 } }};
// Z gate matrix
constexpr dd::Matrix2x2 Zmat = {{{ 1, 0 }, { 0, 0 } }, {{ 0, 0 }, { -1, 0 } }};
// Y gate matrix
constexpr dd::Matrix2x2 Ymat = {{{ 0, 0 }, { 0, -1 } }, {{ 0, 1 }, { 0, 0 } }};
// Hadamard gate matrix
constexpr dd::Matrix2x2 Hmat = {{{ dd::SQRT_2, 0 }, { dd::SQRT_2,  0 }},
                            {{ dd::SQRT_2, 0 }, { -dd::SQRT_2, 0 }}};
// all One. gate matrix (used for test)
constexpr dd::Matrix2x2 Omat = {{{ 1, 0 }, { 1, 0 }},
                                {{ 1, 0 }, { 1, 0 }}};
// number operators
constexpr dd::Matrix2x2 N0mat = {{{ 1, 0 }, { 0, 0 }},
{{ 0, 0 }, { 0, 0 }}};
constexpr dd::Matrix2x2 N1mat = {{{ 0, 0 }, { 0, 0 }},
{{ 0, 0 }, { 1, 0 }}};
class SearchAlgorithms{
    static void BreathFirstTraverese();
};
#endif //DD_PACKAGE_UTIL_H
