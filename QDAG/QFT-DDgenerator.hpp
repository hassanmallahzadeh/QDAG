//
//  QFT-DDgenerator.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-02.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef QFT_DDgenerator_hpp
#define QFT_DDgenerator_hpp
enum PERM_POS//permutation position
{
    NO_PERM = -1,
    BEG_PERM,              // apply first on state
    END_PERN                //apply last on state
};
class StateGenerator{
private:
    dd::Package* dd;
public:
    StateGenerator(dd::Package* dd);
    dd::Edge dd_UniformState(int n);
    dd::Edge dd_Sqrt3State(int n);//if i th digit (before or after fraction point, is divisble by 2 i th bit 0 else 1.
    dd::Edge dd_RandomState(int n);
    dd::Edge dd_BaseState(int n, int i);
};
class GateGenerator{
private:
    dd::Package* dd;
public:
    void lineSet(short* line, int t, int c = -1);
    void lineReset(short* line, int t, int c = -1);
    
    GateGenerator(dd::Package* dd);
    dd::Edge Smatv1(int n, int b1, int b2);
    dd::Edge Smatv2(int n, int b1, int b2);
    void RmatGenerator(dd::Matrix2x2 &m, int k);
    dd::Edge permuteOperator(int n);
};
#endif /* QFT_DDgenerator_hpp */
