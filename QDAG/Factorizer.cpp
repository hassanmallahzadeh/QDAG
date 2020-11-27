//
//  Alice.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#include "Factorizer.hpp"
#include "shorutil.hpp"
#include "ProbabilisticPeriodFinder.hpp"
Factorizer::Factorizer(lli N){
    this->N = N;
}
array<lli, 2> Factorizer::Factors(){
    const int attemptlimit = 100;
    array<lli , 2> factors = {};
    int attemps = 0;
    lli d;//holder for a^(r/2) mod N.
    while(true){
        ++attemps;
        if(attemps > attemptlimit){
            printf("factoring attempts limit (%d) passed, process terminated\n", attemps);
            factors = {-1,-1};
            break;
        }
        lli a = ChooseNumInG_N();
        lli commonfactor = shor::gcd(N,a);
        if(commonfactor != 1){//done, we were lucky!
            factors = {commonfactor, N/commonfactor};
            break;
        }
        lli period = PeriodFinder(a);//TODO: call period finder
        if(period == -1){//period finding limit reached.
            factors = {-1,-1};
            break;
        }
        auto is_odd = [](lli const elem) { return (elem % 2); };
        if(is_odd(period)) continue;
        d = shor::modexp(a,period/2,N);
        if(d + 1 % N != 0) continue;
        lli gcd1 = shor::gcd(d + N - 1/*gcd not defined for 0*/, N);
        lli gcd2 = shor::gcd(d + 1, N);
        factors = {gcd1, gcd2};
        break;
    }
    return factors;
}
lli Factorizer::ChooseNumInG_N(){
    std::random_device device;
    std::mt19937 mt_rand(device());
    std::uniform_int_distribution<lli> dis(2, N - 1);
    return dis(mt_rand);
}
lli Factorizer::PeriodFinder(lli a){
    int trylimit = 100;
    std::pair<lli,lli> p = {-1,-1};
    int counter = 0;
    do{
        if(counter > trylimit){
            printf("period finding attempts limit (%d) passed, process terminated\n", counter);
            return -1;
            break;
        }
        auto* dd = new dd::Package;
        ProbabilisticPeriodFinder pf = ProbabilisticPeriodFinder(N,a,dd);
        p = pf.AttemptReadingMultipleOfInverseOfPeriod();
        if(p.second > 0 && p.first > 0)//if 0 is not measured (input register).
            ++counter;
        delete dd;
    } while(!(p.second>0) || shor::modexp(a,p.second,N) != 1);
    
    return p.second;
}
