//
//  RegisterFactoryInterface.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-21.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "util.h"
#include "RegisterFactory.hpp"
#include "QFT-DDgenerator.hpp"

RegisterFactory::RegisterFactory(lli N, lli a) { 
    this->N = N;
    this->a = a;
}
dd::Edge RegisterFactory::ExponentiateOutReg(){
    dd::Edge e;
    GateGenerator gg = GateGenerator(dd);
    short *line = new short[nt]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < nt; ++i){
        line[i] = -1;
    }
    for(int i = 0; i < ni; ++i){
        gg.lineSet(line, i);
        CtrlMultMod(line, false);
        gg.lineReset(line, i);
        CtrlMultMod(line, true);
    }
    delete[] line;
    return e;
}
void RegisterFactory::CtrlMultMod(const short* line,bool is_rev){
    
}
void RegisterFactory::MakeInitialState(){
    m_state = StateGenerator(dd).dd_BaseState(nt, na + 1);
    GateGenerator gg = GateGenerator(dd);
    short *line = new short[nt]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < nt; ++i){
        line[i] = -1;
    }
    for(int i = 0; i < ni; ++i){
        gg.lineSet(line, i);
        m_state = dd->multiply(dd->makeGateDD(Hmat, nt, line) , m_state);
        gg.lineReset(line, i);
    }
    delete[] line;
};
void RegisterFactory::DetermineRegSizes(){
    for(int i = 0; i < dd::MAXN; ++i){
        long double p = pow(2, i);
        if (ni < 0 && p >= N){
            no = i;
            na = no;
            ni = 2 * i;
            break;
        }
    }
    nt = ni + no + na;
}
