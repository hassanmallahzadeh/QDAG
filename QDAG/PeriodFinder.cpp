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
#include <bitset>
using std::bitset;
//BEG: PeriodFinder start of public methods
/// <#Description#>
/// @param N Number to be factorized, modulo number
/// @param a number to be used as base for finding the period.
PeriodFinder::PeriodFinder(ulli N, ulli a){
    this->N = N;
    this->a = a;
    dd = new dd::Package;
}

//END: PeriodFinder start of public methods

//BEG: PeriodFinder start of private methods

void PeriodFinder::InitializeRegisters(){
    StateGenerator sg = StateGenerator(dd);
    outputRegister = sg.dd_BaseState(no, 0);
    inputRegister = sg.dd_UniformState(ni);
    ulli wreg = a;//work register
    ulli oreg = 1;//num to be put in output register
    vector<dd::ComplexValue> v;
        dd::Edge f = dd::Package::DDone;
           dd::Edge edges[4];
           edges[1] = edges[3] = dd::Package::DDzero;
    ulli workreg = workreg * workreg;//work 'register'
        for (short p = 0; p < no; ++p) {
                   edges[0] = f;
                   edges[2] = f;
               f = dd->makeNonterminal(p, edges);
           }
    int Ntemp = -1;
    while(Ntemp > 0){
        if(Ntemp & 1)
        {
        v.push_back(dd::ComplexValue{1,0});
        oreg *= wreg;
        oreg %= N;
        }
        else
        v.push_back(dd::ComplexValue{0,0});
        wreg = shor::modexp(wreg, 2, N);
        Ntemp>>=1;
    }
     unsigned long remained = ni - v.size();//unsigned long to silence 'precision lost' warning.
    for (int i = 0; i < remained; ++i)
        v.push_back(dd::ComplexValue{0,0});
    inputRegister = sg.dd_CustomState(v, ni);
    v.clear();
    while (oreg > 0) {
        if (oreg & 1)
            v.push_back(dd::ComplexValue{1,0});
        else
            v.push_back(dd::ComplexValue{0,0});
        oreg >>= 1;
    }
    remained = no - v.size();
    for (int i = 0; i < remained; ++i)
         v.push_back(dd::ComplexValue{0,0});
    outputRegister = sg.dd_CustomState(v, no);
}
//END: PeriodFinder start of private methods
