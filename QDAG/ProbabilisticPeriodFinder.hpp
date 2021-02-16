//
//  PeriodFinder.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-08.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef PeriodFinder_hpp
#define PeriodFinder_hpp

#include <stdio.h>
#include "commonheaders.h"
#include "IIC-JKU/DDpackage.h"
#include "QFT-Measurement.hpp"
#include "RegisterFactory.hpp"
class ProbabilisticPeriodFinder{
    lli N;
    lli a;
    dd::Package *dd = nullptr;
    dd::Edge m_state;
    GateGenerator gg;
    void InitializeRegisters();
    vector<bool> m_base2N;//base2 representation of N.
    void HelperRippleAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &);
    void InvHelperRippleAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &);
    void CRippleAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, int t, short *, int, int, dd::Edge &);
    void InvCRippleAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, int t, short *, int, int, dd::Edge &);
     
    void HelperModuloNAdderHalfClassic(const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &, const std::function<int ()> &);
    void CCModuloNAdderHalfClassic(int, int , int, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int , dd::Edge &);
       void InvCCModuloNAdderHalfClassic(int, int , int, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int , dd::Edge &);
    void CCRippleAdderHalfClassic(int, int, int, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &);
      void InvCCRippleAdderHalfClassic(int, int, int, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, const std::function<int (int)> &, short *, int, dd::Edge &);
    void HelperCCRippleHalfClassic(std::function<void (int, int, int, int)> &, std::function<void (int, int, int, int)> &, std::function<void (int, int, int)> &, short *, int, int, dd::Edge &, int);
    void CMultiplierModuloNHalfClassic(const std::function<int (int)> &a0Nbase2, const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, lli a0cnum, short *line, int mcindex, int nt, dd::Edge &state, const std::function<int ()> &tindex, const std::function<int (int)> &xindice);
   void InvCMultiplierModuloNHalfClassic(const std::function<int (int)> &a0Nbase2, const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, lli a0cnum, short *line, int mcindex, int nt, dd::Edge &state, const std::function<int ()> &tindex, const std::function<int (int)> &xindice);
public:
    
    static const int posControl;
    int n = -1;//output register size
    std::pair<lli,lli> AttemptReadingMultipleOfInverseOfPeriod();
    ProbabilisticPeriodFinder(lli, lli, dd::Package * = nullptr);
    ~ProbabilisticPeriodFinder();
};
#endif /* PeriodFinder_hpp */
