//
//  PeriodFinder.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#include "PeriodFinder.hpp"
#include "QFT-DDgenerator.hpp"
#include "shorutil.hpp"
#include "QFT-Measurement.hpp"
#include <bitset>
using std::bitset;
//BEG: PeriodFinder start of public methods
/// Period Finder Constructor.
/// @param N Number to be factorized, modulo number
/// @param a number to be used as base for finding the period.
PeriodFinder::PeriodFinder(lli N, lli a){
    this->N = N;
    this->a = a;
    dd = new dd::Package;
    RegisterFactory rf(N, a, dd);
}

//END: PeriodFinder start of public methods

//BEG: PeriodFinder start of private methods

void PeriodFinder::InitializeRegisters(){
   vector<bool> base2N = shor::base2rep(N);//base 2 representation of N
    no = static_cast<int>(base2N.size());//num qubits needed to represent N.
    ni = 2 * no;//input register size.
    vector<int> numq;
    for(int i = 0; i < ni; ++i){
        numq.push_back(2);//2 means all qubits in numq will be (|0> + |1>)/2
    }
    state = rf.ExponentiatorModuloN(numq);
}

void PeriodFinder::MeasureOutputReg(){
    std::function<int (int)> outputindice = rf.OutPutRegIndice();
    Measurement mm(dd);
    std::random_device device;
    std::mt19937 mt_rand(device);
    for(int i = 0; i < no; ++i){
        mm.Measure(state, rf.nt, outputindice(i), mt_rand);
    };
}
dd::Edge PeriodFinder::DebugPeriodFinder(){
    return state;
}
PeriodFinder::~PeriodFinder(){
    delete dd;
}
//END: PeriodFinder start of private methods
