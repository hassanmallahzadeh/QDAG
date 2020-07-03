//
//  QFT.hpp
//  cmakelearn
//
//  Created by Hassan Mallahzadeh on 2020-06-09.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef QFT_hpp
#define QFT_hpp

#include <stdio.h>
enum PERM_POS//permutation position
{
    NO_PERM = -1,
    BEG_PERM,              // apply first on state
    END_PERN                //apply last on state
};
dd::Edge Smatv1(dd::Package *dd, int n, int b1, int b2);
dd::Edge Smatv2(dd::Package *dd, int n, int b1, int b2);
void RmatGenerator(dd::Matrix2x2 &m, int k);
dd::Edge permuteOperator(dd::Package *dd, int n);
void lineSet(short* line, int t, int c = -1);
void lineReset(short* line, int t, int c = -1);
dd::Edge dd_Inverse(dd::Package *dd, int n);
void recordDDsize(const dd::Edge& e, std::ofstream& datafile, int v, int g = -1);
dd::Edge dd_QFTV1(dd::Package *dd, int n, PERM_POS perm);
dd::Edge dd_QFTV2(dd::Package *dd, int n, PERM_POS perm);
#endif /* QFT_hpp */
