//
//  Bob.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#include "Bob.hpp"

void Bob::ForcePandQ(ulli p, ulli q) {
    //assert(make sure they are primes)//TODO: nice to make sure primes are passed as p and q.
    this->p = p;
    this->q = q;
}
ulli Bob::ReturnN(){
    return p * q;
}
array<ulli, 2> Bob::randPrimes(){
    assert(0);//TODO: nice to have ability to generate random primes. Care must be taken not to choose 'unsafe' N = pq as N = 3*5
    array<ulli, 2> ret;
    return ret;
}
