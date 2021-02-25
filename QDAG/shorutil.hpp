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
vector<bool> base2rep(lli N, int n = -1);
std::pair<lli,lli> contfrac(lli, int, int);
lli base2to10(vector<int>,bool);
lli pow(lli,lli);
lli bfpf(lli,lli);
}
#endif /* shorutil_hpp */
