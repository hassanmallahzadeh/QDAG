//
//  RegisterFactory.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-21.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "util.h"
#include "RegisterFactory.hpp"
#include "QFT-DDgenerator.hpp"

RegisterFactory::RegisterFactory(lli N, lli a, dd::Package *dd) : gg(dd) {
    this->N = N;
    this->a = a;
    this->dd = dd;
    line = new short[nt]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < nt; ++i){
        line[i] = -1;
    }
}
/// CExponentiation the input register, put in output register.
dd::Edge RegisterFactory::CExponentiation(){
    dd::Edge e;
    short *line = new short[nt]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < nt; ++i){
        line[i] = -1;
    }
    for(int i = 0; i < ni; ++i){
        gg.lineSet(line, i);
        CtrlMultMod(line, false);
        gg.permuteOperatorOnState(nt, ni, no, na, m_state);
        CtrlMultMod(line, true);
        gg.lineReset(line, i);
    }
    delete[] line;
    return e;
}
void RegisterFactory::CtrlMultMod(const short* line,bool is_rev){
    
}

void RegisterFactory::MakeInitialState(){
    m_state = StateGenerator(dd).dd_BaseState(nt, na + 1);
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
void RegisterFactory::Carry(int c0, int a0, int b0, int c1){
//    m_state = dd->multiply(, m_state);
}
void RegisterFactory::ToffoliSet(int c0 = -1, int c1 = -1, int t = -1){
    if(c0 > 0)
        line[c0] = 0;
}
RegisterFactory::~RegisterFactory(){
    delete[] line;
}
