//
//  QFT.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-06-09.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
///

#ifndef QFT_hpp
#define QFT_hpp
#include "commonheaders.h"
#include "IIC-JKU/DDpackage.h"
#include "QFT-DDgenerator.hpp"

#include <stdio.h>
class QFT{
private:
    dd::Package *dd = nullptr;
    enum ORD_C_T{
        REG_C_T,//figure 3.1 Mermin. regular control/target setting.
        INV_C_T,//figure 3.2 Mermin. inverse control/target setting.
    };
    ORD_C_T m_ord = REG_C_T;//change to invert control/target setting
public:
    static const int posControl;
    QFT(dd::Package *dd);
    dd::Edge dd_QFTV1(int n,  PERM_POS perm);
    dd::Edge dd_QFTV2(int n, dd::Edge state, PERM_POS perm, vector<int> v, bool = false);
    dd::Edge dd_QFTV3(int n,  PERM_POS perm);
    dd::Edge dd_QFTV4(int n, dd::Edge state, PERM_POS perm);
    void dd_QFTV5(int n, dd::Edge &state, vector<int> &indice, bool inv);//TODO: pass state by reference and cahnge directly.
    void dd_QFTV5Reverse(int n, dd::Edge &state, vector<int> indice, bool inv);
    dd::Edge dd_QFTGNV1(int n, dd::Edge state, PERM_POS perm, engine& unrg);//Griffiths–Niu
    dd::Edge dd_QFTGNV2(int n, dd::Edge state, PERM_POS perm, engine& urng);//Griffiths–Niu,
    vector<int> dd_QFTGNV1(int n, dd::Edge &state, PERM_POS perm, engine& unrg, vector<int> v);//Griffiths–Niu
};
#endif /* QFT_hpp */
