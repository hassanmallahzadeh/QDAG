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
//parameters: variable roles in gate(to be set, number of variables, target index, control index(-1:
//none)
void lineSet(short* line, int t, int c = -1);
void lineReset(short* line, int t, int c = -1);
dd::Edge dd_Inverse(dd::Package *dd, int n);
void recordDDsize(const dd::Edge& e, std::ofstream& datafile, int v, int g = -1);
#endif /* QFT_hpp */
