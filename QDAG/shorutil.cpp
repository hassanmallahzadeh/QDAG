//
//  shorutil.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#include "shorutil.hpp"
namespace shor{
/// greatest common divisor, Euclidean algorithm.
/// @param f first number
/// @param c second number
lli gcd(lli f, lli c){
    assert( f > 0 & c > 0);
    if(f == c){
        return f;
    }
    if(f < c){
        std::swap(f,c);
    }
    lli ftemp = -1;
    while(c > 0){
        ftemp = f;
        f = c;
     //   c = ftemp - floor(ftemp/c) * c;
        c = ftemp % c;
    }
    return f;
}



/// This function takes three integers, x, a, and n, and returns x^a mod n.  This algorithm is known as the "Russian peasant method," I believe, and avoids overflow by never calculating x^a directly.
/// @param x base
/// @param a exponent
/// @param n mod
lli modexp( lli x,  lli a,  lli n) {
   lli value = 1;
   lli tmp = x % n;
  while (a > 0) {
    if (a & 1) {
      value = (value * tmp) % n;
    }
    tmp = tmp * tmp % n;
    a>>=1;
  }
  return value;
}
}
