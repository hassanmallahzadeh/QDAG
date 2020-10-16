//
//  Bob.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef Bob_hpp
#define Bob_hpp
#include "commonheaders.h"
#include <stdio.h>
class Bob{
    lli p = -1;
    lli q = -1;
    array<lli, 2> randPrimes();
public:
    void ForcePandQ(lli p, lli q);
    lli ReturnN();
};
#endif /* Bob_hpp */
