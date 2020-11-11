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
#include "QFT-Measurement.hpp"
#include "RegisterFactory.hpp"
#include <bit>
class PeriodFinder{
    lli N;
    lli a;
    int ni = -1;//input register size
    int no = -1;//output register size
    dd::Package *dd = nullptr;
    dd::Edge state;
    RegisterFactory* p_rf;
    void InitializeRegisters();
    void MeasureOutputReg();
    lli ApplyQFT();
    lli AttemptFindingPeriod();
public:
    std::pair<lli,lli> DebugPeriodFinder();
    PeriodFinder(lli, lli, dd::Package *);
    ~PeriodFinder();
};
#endif /* PeriodFinder_hpp */
