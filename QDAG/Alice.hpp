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
    ulli N;
    ulli ChooseNumInG_N();
    engine unrg;
    const int attemptlimit = 100;
public:
    Alice(ulli N, engine& unrg);
    array<ulli, 2> Factors();
};
#endif /* Alice_hpp */
