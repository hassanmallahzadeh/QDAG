//
//  PeriodFinder.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef PeriodFinder_hpp
#define PeriodFinder_hpp

#include <stdio.h>
#include "commonheaders.h"
#include "IIC-JKU/DDpackage.h"
#include <bit>
class PeriodFinder{
    ulli N;
    ulli a;
    int ni = -1;//input register size
    int no = -1;//output register size
    dd::Package *dd = nullptr;
    dd::Edge inputRegister = {nullptr, {nullptr, nullptr}};
    dd::Edge outputRegister = {nullptr, {nullptr, nullptr}};
    void InitializeRegisters();
public:
    PeriodFinder(ulli N, ulli b);
    ulli FoundPeriod();
    bool IsNumberPow2();
};
#endif /* PeriodFinder_hpp */
