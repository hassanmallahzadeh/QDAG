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
class Alice{
    lli N;
    lli ChooseNumInG_N();
    engine unrg;
    const int attemptlimit = 100;
public:
    Alice(lli N, engine& unrg);
    array<lli, 2> Factors();
};
#endif /* Alice_hpp */
