//
//  QFT-DDgenerator.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-02.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef QFT_DDgenerator_hpp
#define QFT_DDgenerator_hpp
#include "commonheaders.h"
#include <random>
using CN = dd::ComplexNumbers;
enum PERM_POS//permutation position
{
    NO_PERM = -1,
    BEG_PERM,              // apply first on state: original variable order
    END_PERM                //apply last on state: inverse variable order (?)
};
class StateGenerator{
private:
    dd::Package* dd;
public:
    StateGenerator(dd::Package* dd);
    dd::Edge dd_UniformState(int n);
    dd::Edge dd_Sqrt3State(int n);//if i th digit (before or after fraction point, is divisble by 2 i th bit 0 else 1.
    dd::Edge dd_RandomState(int n, int seed);
    dd::Edge dd_BaseState(int n, lli i);
    dd::Edge dd_CustomState(vector<dd::ComplexValue> v, int n);
};
class GateGenerator{
private:
    dd::Package* dd;
public:
    void lineSet(short* line, int t, int c0 = -1, int c1 = -1, bool c0p = true, bool c1p = true);
    void lineSet(short* line, int t, const map<int,bool> &m);
    void lineReset(short* line, int t, int c0 = -1, int c1 = -1);
    void lineReset(short* line, int t, const map<int,bool> &m);
    static void lineClear(short* line, int n);
    GateGenerator(dd::Package* dd);
    dd::Edge Smatv1(int n, int b1, int b2);
    dd::Edge Smatv2(int n, int b1, int b2);
    void RmatGenerator(dd::Matrix2x2 &m, int k);
    void RInvmatGenerator(dd::Matrix2x2 &m, int k);
    dd::Edge permuteOperator(int n);
    dd::Edge permuteOperatorOnState(int n, dd::Edge state);
    dd::Edge swapRegistersOnState(int nt, vector<int> v1, vector<int> v2, dd::Edge state);
    dd::Edge ToffoliGenOrApply(short* line, int, int, int, int, dd::Edge* state = nullptr ,bool = true, bool = true);
    dd::Edge CNotGenOrApply(short* line, int, int, int, dd::Edge* state = nullptr, bool = true);
    dd::Edge NotGenOrApply(short* line, int, int, dd::Edge* state = nullptr);
    dd::Edge CKNotGenOrApply(short* line, int t, const map<int,bool> &m, int nt, dd::Edge* state = nullptr);
    dd::Edge HadGenOrApply(short* line, int, int, dd::Edge* state = nullptr);
    dd::Edge RmatGenOrApply(short* line, int, int, int, dd::Edge* state = nullptr);
    dd::Edge CKRmatGenOrApply(short* line, int , int , const map<int,bool> & , int , dd::Edge* = nullptr, bool isinverse = false);
};
#endif /* QFT_DDgenerator_hpp */
