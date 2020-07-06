//
//  QFT.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-06-09.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
///

#ifndef QFT_hpp
#define QFT_hpp
#include "IIC-JKU/DDpackage.h"
#include "QFT-DDgenerator.hpp"
#include <stdio.h>
class QFT{
private:
    dd::Package *dd = nullptr;
public:
    QFT(dd::Package *dd);
    dd::Edge dd_QFTV1(int n,  PERM_POS perm);
    dd::Edge dd_QFTV2(int n, dd::Edge state, PERM_POS perm);
};
#endif /* QFT_hpp */
