//
//  RegisterFactoryInterface.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-21.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef RegisterFactoryInterface_hpp
#define RegisterFactoryInterface_hpp
#include "commonheaders.h"
#include "DDpackage.h"
#include "QFT-DDgenerator.hpp"
#include <stdio.h>
class RegisterFactory{
private:
    lli N = -1;
    lli a = -1;
    int ni = -1;//input register size
    int no = -1;//output register size
    int na = -1;//auxiliary register size
    dd::Edge m_state = dd::Package::DDnull;
    dd::Package *dd = nullptr;
    int nt = -1;//total number of qubits (input+output+auxiliary)
    short* line = nullptr;//holds 'line'for setting basic gates.
    GateGenerator gg;
    
    void DetermineRegSizes();
    void MakeInitialState();
    void CtrlMultMod(const short* = nullptr, bool = false);
    void Carry(int c0,int a0, int b0, int c1);
    void ToffoliSet(int , int, int);
    void ToffoliReset(int, int, int);
public:
    RegisterFactory(lli N, lli a, dd::Package *dd);
    ~RegisterFactory();
    dd::Edge CExponentiation();
};
#endif /* RegisterFactoryInterface_hpp */
