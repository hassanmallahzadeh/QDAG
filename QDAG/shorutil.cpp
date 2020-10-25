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
    lli ftemp = 0;
    while(c > 0){
        ftemp = f;
        f = c;
        //   c = ftemp - floor(ftemp/c) * c;
        c = ftemp % c;
    }
    return f;
}
// Returns modulo inverse of a with respect
// to m using extended Euclid Algorithm
// Assumption: a and m are coprimes, i.e.,
// gcd(a, m) = 1//https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/
lli modInverse(lli a, lli m)
{
    assert(gcd(a,m) == 1);//they must be co-prime.
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
        static_cast<void>(m = a % m), a = t;
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

int bddNumVar(dd::Edge edge, bool isVector){
    list<dd::NodePtr> l;
    assert(edge.p);// make sure not null
    dd::NodePtr curnode = edge.p;
    l.push_back(curnode);
    int max = 0;
    while(!l.empty()){// breath first search
        curnode = l.front();
        if(curnode->v > max){
            max = curnode->v;
        }
        l.pop_front();
        for(int i = 0; i < dd::NEDGE; ++i){
            if(isVector && i % dd::RADIX != 0) continue;
            if(curnode->e[i].p == dd::Package::terminalNode){// base case
                continue;
            }
            else{
                l.push_back(curnode->e[i].p);
            }
        }
    }
    return max + 1;
}
//returns base 2 representation. returns in reverse by default.
/// Constructs base2 representation of N, n digits. If n > log_2(N) + 1, pad with zeros
/// @param N Number for which base2 rep is to be found
/// @param n number of digits (bits) available. most be n > log_2(N).
vector<bool> base2rep(lli N,int n){
    vector<bool> v;
    int dc = 0;//digit counter
    while(N > 0){
        v.push_back(N & 1);
        N >>= 1;
        ++dc;
    }
    if(n > 0){
    assert(n >= dc);
    for(int i = 0; i < n - dc; ++i){
        v.push_back(0);
    }
    }
    return v;
}

}
