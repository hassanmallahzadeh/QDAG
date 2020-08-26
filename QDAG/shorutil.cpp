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
// From geeksforgeeks.com
// Iterative C++ program to find modular
// inverse using extended Euclid algorithm
// Returns modulo inverse of a with respect
// to m using extended Euclid Algorithm
// Assumption: a and m are coprimes, i.e.,
// gcd(a, m) = 1
/// modular multtiplication inverse of 'a' mod m.
/// @param a <#a description#>
/// @param m, N =pq in our problem
lli modInverse(lli a, lli m)
{
    if(m < a){
        std::swap(m,a);
    }
    lli m0 = m;
    lli y = 0, x = 1;
  
    if (m == 1)
      return 0;
  
    while (a > 1)
    {
        // q is quotient
        lli q = a / m;
        lli t = m;
  
        // m is remainder now, process same as
        // Euclid's algo
        m = a % m;
        a = t;
        t = y;
  
        // Update y and x
        y = x - q * y;
        x = t;
    }
  
    // Make x positive
    if (x < 0)
       x += m0;
  
    return x;
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
