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
#include "IIC-JKU/DDpackage.h"
#include "QFT-DDgenerator.hpp"
#include <stdio.h>
class RegisterFactory{
private:
    ulli N = 0;
    ulli a = 0;
    int n = -1;
    dd::Edge state = dd::Package::DDnull;
    dd::Package *dd = nullptr;
    vector<bool> base2N;//base2 representation of N.
    short* line = nullptr;//holds 'line'for setting basic gates.
    GateGenerator gg;
public:
    RegisterFactory(ulli N, ulli a, dd::Package *dd);
    dd::Edge RippleAdderDebug(vector<int> num1, vector<int> num2);
    std::function<void (int, int)> StateInitializer(short *line, int nt, dd::Edge &state);
    dd::Edge RippleAdderHalfClassicDebug(ulli cl, vector<int> num);
    void HelperRippleAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &);
    void InvHelperRippleAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &);
    void CRippleAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, int t, short *, int, int, dd::Edge &);
    void InvCRippleAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, int t, short *, int, int, dd::Edge &);
     
    void HelperModuloNAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &, const std::function<int ()> &);
    void CCModuloNAdderHalfClassic(int, int , int, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int , dd::Edge &);
      void InvCCModuloNAdderHalfClassic(int, int , int, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int , dd::Edge &);
    dd::Edge ModuloNAdderHalfClassicDebug(ulli cnum, vector<int> qnum);
    void CMultiplierModuloNHalfClassic(const std::function<int (int)> &a0Nbase2, const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, ulli a0cnum, short *line, int mcindex, int nt, dd::Edge &state, const std::function<int ()> &tindex, const std::function<int (int)> &xindice);
    void InvCMultiplierModuloNHalfClassic(const std::function<int (int)> &a0Nbase2, const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, ulli a0cnum, short *line, int mcindex, int nt, dd::Edge &state, const std::function<int ()> &tindex, const std::function<int (int)> &xindice);
       
    dd::Edge CMultiplierModuloNDebug(ulli cnum, vector<int> qnum, int mcv);
    dd::Edge ExponentiatorModuloNDebug(ulli, vector<int>);
    void HelperCCRippleHalfClassic(std::function<void (int, int, int, int)> &, std::function<void (int, int, int, int)> &, std::function<void (int, int, int)> &, short *, int, int, dd::Edge &, int);
    
    void CCRippleAdderHalfClassic(int, int, int, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &);
      void InvCCRippleAdderHalfClassic(int, int, int, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &);
    
    ~RegisterFactory();
    dd::Edge CExponentiation();
};
#endif /* RegisterFactoryInterface_hpp */
