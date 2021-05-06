//
//  SB-ProbabilisticPeriodFinder.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2021-03-18.
//  Copyright © 2021 Hassan Mallahzadeh. All rights reserved.
//

#ifndef SB_ProbabilisticPeriodFinder_hpp
#define SB_ProbabilisticPeriodFinder_hpp
#include "commonheaders.h"
#include "IIC-JKU/DDpackage.h"
#include "QFT-Measurement.hpp"
#include "RegisterFactory.hpp"
#include "QFT.hpp"
#include <stdio.h>
class SB_PPF{
    lli N;
    lli a;
    int nt;
    dd::Package *dd = nullptr;
    dd::Edge m_state;
    GateGenerator gg;
    QFT m_qft;
    void phi_Adder(const std::function<int (int)> &acnum, const std::function<int (int)> &bindice, const map<int,bool> &c, short *line, dd::Edge &state);
    void phi_Subtractor(const std::function<int (int)> &acnum, const std::function<int (int)> &bindice, const map<int,bool> &c, short *line, dd::Edge &state);
    void phi_AdderModN(const std::function<int (int)> &acnum, const std::function<int (int)> &bindice, const std::function<int ()> &tindex, const map<int,bool> &c, short *line, dd::Edge &state);
    void phi_SubtractorModN(const std::function<int (int)> &acnum, const std::function<int (int)> &bindice, const std::function<int ()> &tindex, const map<int,bool> &c, short *line, dd::Edge &state);
    void phi_CMultiplier(const lli &acnum, const std::function<int (int)> &bindice, const std::function<int (int)> &xindice, const std::function<int ()> &tindex, map<int,bool> c, short *line, dd::Edge &state);
    void phi_CDivider(const lli &acnum, const std::function<int (int)> &bindice, const std::function<int (int)> &xindice, const std::function<int ()> &tindex, map<int,bool> c, short *line, dd::Edge &state);
    void phi_CUa(const lli &acnum, const std::function<int (int)> &bindice, const std::function<int (int)> &xindice, const std::function<int ()> &tindex, map<int,bool> c, short *line, dd::Edge &state);
    std::function<int (int)> m_base2N;//wrapper to read N classic bits.
//    void out_Circuit(const std::function<int (int)> &bindice, const std::function<int (int)> &xindice, const std::function<int ()> &tindex, , const std::function<int ()> &mindex, map<int,bool> c, short *line, dd::Edge &state);
    dd::Matrix2x2 Rmat;// rotation gate
public:
    const int posControl;//apply rotation gates on stage i if qubit measrement on stage i had result posControl. refer to Mermin, fig 3.3. Would be 0, 1 or 2 for base 3, etc.
    int n = -1;//output register size
    std::pair<lli,lli> AttemptReadingMultipleOfInverseOfPeriod();
    
    SB_PPF(lli, lli, dd::Package * = nullptr);
   // ~SB_PPF();
};
#endif /* SB_ProbabilisticPeriodFinder_hpp */
