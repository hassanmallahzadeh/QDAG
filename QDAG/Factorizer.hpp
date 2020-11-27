//
//  Alice.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef Alice_hpp
#define Alice_hpp
#include "commonheaders.h"
#include <stdio.h>
class Factorizer{
    lli N;
    lli ChooseNumInG_N();
public:
    Factorizer(lli N);
    array<lli, 2> Factors();
    lli PeriodFinder(lli a);
};
#endif /* Alice_hpp */
