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
PeriodFinder::PeriodFinder(lli N, lli a, dd::Package *dd){
    this->N = N;
    this->a = a;
    this->dd = dd;
    p_rf = new RegisterFactory(N, a, dd);
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
    state = p_rf->ExponentiatorModuloN(numq);
}

void PeriodFinder::MeasureOutputReg(){
    std::function<int (int)> outputindice = p_rf->OutPutRegIndice();
    
    std::random_device device;
    std::mt19937 mt_rand(device());

    dd->export2Dot(state, "before.dot");
    for(int i = 0; i < no; ++i){
        //TODO: measurement mechanism is 'reset' for every qubit. Look into reusing former calculations.
        Measurement mm(dd);
        int res = mm.Measure(state, p_rf->nt, outputindice(i), mt_rand);
        dd->garbageCollect();
        std::cout<<"qubit "<<outputindice(i)<<" result: "<<res<<std::endl;
        dd->export2Dot(state, "iter" + std::to_string(i) + ".dot");
    };
}
dd::Edge PeriodFinder::DebugPeriodFinder(){
    InitializeRegisters();
    MeasureOutputReg();
    return state;
}
PeriodFinder::~PeriodFinder(){
    delete p_rf;
}
//END: PeriodFinder start of private methods
