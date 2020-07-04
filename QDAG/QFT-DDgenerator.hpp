//
//  QFT-DDgenerator.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-02.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef QFT_DDgenerator_hpp
#define QFT_DDgenerator_hpp
#include <util.h>
class StateGenerator{
private:
    dd::Package* dd;
public:
    StateGenerator();
    dd::Edge dd_UniformState(int n);
    dd::Edge dd_Sqrt3State(int n);//if i th digit (before or after fraction point, is divisble by 2 i th bit 0 else 1.
    dd::Edge dd_RandomState(int n);
    dd::Edge dd_BaseState(int n, int i);
};
#endif /* QFT_DDgenerator_hpp */
