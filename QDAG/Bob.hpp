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
    ulli p = -1;
    ulli q = -1;
    array<ulli, 2> randPrimes();
public:
    void ForcePandQ(ulli p, ulli q);
    ulli ReturnN();
};
#endif /* Bob_hpp */
