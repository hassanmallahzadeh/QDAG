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
lli gcd(lli f, lli c);
lli modexp( lli x,  lli a,  lli n);
lli modInverse(lli a, lli m);
int bddNumVar(dd::Edge root, bool isVector);
vector<bool> base2rep(lli N, bool inversed = true);
}
#endif /* shorutil_hpp */
