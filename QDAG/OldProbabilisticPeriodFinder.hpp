//
//  PeriodFinder.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef OldPeriodFinder_hpp
#define OldPeriodFinder_hpp

#include <stdio.h>
#include "commonheaders.h"
#include "IIC-JKU/DDpackage.h"
#include "QFT-Measurement.hpp"
#include "RegisterFactory.hpp"
class OldProbabilisticPeriodFinder{
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
    OldProbabilisticPeriodFinder(lli, lli, dd::Package * = nullptr);
    ~OldProbabilisticPeriodFinder();
};
#endif /* PeriodFinder_hpp */

