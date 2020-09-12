//
//  Alice.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#include "Alice.hpp"
#include "shorutil.hpp"
//BEG: Alice Public Methods
Alice::Alice(ulli N, engine &unrg){
    this->N = N;
    this->unrg = unrg;
}
array<ulli, 2> Alice::Factors(){
    array<ulli , 2> factors = {};
    int attemps = 0;
    ulli d;//holder for a^(r/2) mod N.
    while(true){
        ++attemps;
        if(attemps > attemptlimit){
            std::cout<<"Max factoring attempts reached\n";
            assert(0);
        }
        ulli a = ChooseNumInG_N();
        ulli commonfactor = shor::gcd(N,a);
        if(commonfactor != 1){//done, we were lucky!
            factors = {commonfactor, N/commonfactor};
            break;
        }
        ulli period = -1;//TODO: call period finder
        auto is_odd = [](ulli const elem) { return (elem % 2); };
        if(is_odd(period)){
            continue;
        }
        else{
            d = shor::modexp(a,period/2,N);
            if(d + 1 % N != 0)
                continue;
            ulli gcd1 = shor::gcd(d + N - 1/*gcd not defined for 0*/, N);
            ulli gcd2 = shor::gcd(d + 1, N);
            factors = {gcd1, gcd2};
            break;
        }
    }
    return factors;
}
//END: Alice Public Methods

//BEG: Alice Private Methods
ulli Alice::ChooseNumInG_N(){
    std::uniform_int_distribution<ulli> dis(N - 1);
    return dis(unrg);
}
//END: Alice Private Methods
