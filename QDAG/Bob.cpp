//
//  Bob.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#include "Bob.hpp"

void Bob::ForcePandQ(lli p, lli q) {
    //assert(make sure they are primes)//TODO: nice to make sure primes are passed as p and q.
    this->p = p;
    this->q = q;
}
lli Bob::ReturnN(){
    return p * q;
}
array<lli, 2> Bob::randPrimes(){
    assert(0);//TODO: nice to have ability to generate random primes. Care must be taken not to choose 'unsafe' N = pq as N = 3*5
    array<lli, 2> ret;
    return ret;
}
