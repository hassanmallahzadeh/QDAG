//
//  shorutil.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef shorutil_hpp
#define shorutil_hpp
#include "commonheaders.h"
#include "IIC-JKU/DDpackage.h"
#include <stdio.h>
namespace shor{
ulli gcd(ulli f, ulli c);
ulli modexp( ulli x,  ulli a,  ulli n);
ulli modInverse(ulli a, ulli m);
int bddNumVar(dd::Edge root, bool isVector);
vector<bool> base2rep(ulli N, int n = -1);
}
#endif /* shorutil_hpp */
