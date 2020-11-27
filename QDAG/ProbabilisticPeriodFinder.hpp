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
class ProbabilisticPeriodFinder{
    lli N;
    lli a;
    dd::Package *dd = nullptr;
    dd::Edge state;
    RegisterFactory* p_rf;
    void InitializeRegisters();
    lli MeasureOutputReg();
    lli ApplyQFTMeasureInputRegister();
  
public:
    int ni = -1;//input register size
    int no = -1;//output register size
    std::pair<lli,lli> AttemptReadingMultipleOfInverseOfPeriod();
    ProbabilisticPeriodFinder(lli, lli, dd::Package * = nullptr);
    ~ProbabilisticPeriodFinder();
};
#endif /* PeriodFinder_hpp */
