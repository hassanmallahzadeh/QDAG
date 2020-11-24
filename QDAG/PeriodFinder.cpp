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
#include "QFT.hpp"
//#define DEMONSTRATE
//BEG: PeriodFinder start of public methods
/// Period Finder Constructor.
/// @param N Number to be factorized, modulo number
/// @param a number to be used as base for finding the period.
PeriodFinder::PeriodFinder(lli N, lli a, dd::Package *dd){
    this->N = N;
    this->a = a;
    dd ? this->dd = dd : dd = new dd::Package;
    vector<bool> base2N = shor::base2rep(N);//base 2 representation of N
     no = static_cast<int>(base2N.size());//num qubits needed to represent N.
     base2N = shor::base2rep(N*N);//base 2 representation of N^2
     ni = static_cast<int>(base2N.size());//num qubits needed to represent N^2.
    p_rf = new RegisterFactory(N, a, ni, no, dd);
}

//END: PeriodFinder start of public methods

//BEG: PeriodFinder start of private methods

void PeriodFinder::InitializeRegisters(){

    vector<int> numq;
    for(int i = 0; i < ni; ++i){
        numq.push_back(2);//2 means all qubits in numq will be (|0> + |1>)/2
    }
    state = p_rf->ExponentiatorModuloN(numq);
}

lli PeriodFinder::MeasureOutputReg(){
    bool report = false;
    std::function<int (int)> outputindice = p_rf->OutPutRegIndice();
    
    std::random_device device;
    std::mt19937 mt_rand(device());
    
    if(report)
    std::cout<<"measurement on output register after exponentiator:\n";
    vector<int> vres;
    for(int i = 0; i < no; ++i){
        //TODO: measurement mechanism is 'reset' for every qubit. Look into reusing former calculations.
        Measurement mm(dd);
        int res = mm.Measure(state, p_rf->nt, outputindice(i), mt_rand);
        vres.push_back(res);
if(report)
        std::cout<<"qubit "<<outputindice(i)<<" result: "<<res<<std::endl;
    };
    return shor::base2to10(vres, false);
}
/// Apply QFT to inpput register, perform emasurement with GN scheme, return resulted number in base 10
///@return number outcome of measurement on input register.
lli PeriodFinder::ApplyQFT(){
    std::function<int (int)> inputindice;
    inputindice = p_rf->InputRegIndice();
    QFT *p_qft;//I do on heap since this is a bdd...
    std::mt19937 mt_rand((std::random_device())());
    p_qft = new QFT(dd);
    vector<int> indices;
    for(int i = 0; i < ni; ++i){
        indices.push_back(inputindice(i));
    }
    std::reverse(indices.begin(), indices.end());//TAKE OUT DEBUG
    vector<int> qftgnmo;//qft griffiths niu measurement outcomes.
    qftgnmo = p_qft->dd_QFTGNV1(p_rf->nt, state, PERM_POS::NO_PERM, mt_rand, indices);
    assert(qftgnmo.size() == ni);
//    std::reverse(qftgnmo.begin(), qftgnmo.end());//QFT has reversed output.
#ifdef DEMONSTRATE
    std::cout<<"measurement on input register during qft:\n";
    for(int i = 0; i < indices.size(); ++i){
    std::cout<<"qubit "<<indices[i]<<" result: "<<qftgnmo[i]<<std::endl;
    }
    std::cout<<"in base 10:\n";
    std::cout << y << std::endl;
#endif
    lli y = shor::base2to10(qftgnmo, false);
    delete p_qft;
    return y;
}
std::pair<lli,lli> PeriodFinder::DebugPeriodFinder(){
    InitializeRegisters();
    MeasureOutputReg();
    lli inregout = ApplyQFT();
    std::pair<lli,lli> p = shor::contfrac(inregout, ni, no);
    return p;
}
std::pair<lli,lli> PeriodFinder::DebugMeasureInputRegAfterQFT(){
    InitializeRegisters();
    lli ores = MeasureOutputReg();
    lli ires = ApplyQFT();
    return std::pair<lli,lli>(ires,ores);;
}
std::pair<lli,lli> PeriodFinder::DebugMeasureInputRegNoQFT(){
    InitializeRegisters();
    lli ores = MeasureOutputReg();
    auto inputindice = p_rf->InputRegIndice();
    std::mt19937 mt_rand((std::random_device())());
    vector<int> vres;
    for(int i = 0; i < ni; ++i){
        //TODO: measurement mechanism is 'reset' for every qubit. Look into reusing former calculations.
        Measurement mm(dd);
        int res = mm.Measure(state, p_rf->nt, inputindice(i), mt_rand);
        vres.push_back(res);
    };
    lli ires = shor::base2to10(vres, false);
    return std::pair<lli,lli>(ires,ores);
}
PeriodFinder::~PeriodFinder(){
    delete p_rf;
}
//END: PeriodFinder start of private methods
