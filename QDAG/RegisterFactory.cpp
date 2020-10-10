//
//  RegisterFactory.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-21.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved. Implemented based on PRA S1050-2947(96)05707-1, refed to as the 'paper' in code.
//
#include "util.h"
#include "RegisterFactory.hpp"
#include "QFT-DDgenerator.hpp"
#include "shorutil.hpp"
#include <functional>

/// RegisterFactory constructor
/// @param N number to be factored
/// @param a random number smaller than N
/// @param dd pointer to Package class
RegisterFactory::RegisterFactory(ulli N, ulli a, dd::Package *dd) : gg(dd) {
    this->dd = dd;
    if(N > 0){//0 blocks for debugging of ripple carry
        this->N = N;//num to be factorized
        this->a = a;//randomly chosen number smaller than N
        base2N = shor::base2rep(N);//base 2 representation of N
        n = static_cast<int>(base2N.size());//num qubits needed to represent N.
    }
}
/// Inverse of Controlled Multiplication Mod N.
/// Used to see if ripple adder works as expected. Based on figure 2 of the paper.
///@param num1 first number in addition. size is base 2 representation. for each element, 0 means variable = 0, 1 means varibale = 1, 2 means (|0> + |1>)√2, 3 means |0> - |1>)√2
///@param num2 second number in addition, description like num1
dd::Edge RegisterFactory::RippleAdderDebug(vector<int> num1, vector<int> num2){
    assert(num1.size() == num2.size());
    assert (n == static_cast<int> (num1.size()));
    
    int nt = n * 3 + 1;//first number, second number, carry and overflow carry.
    short* line = new short[nt];
    gg.lineClear(line, nt);
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    std::function<void (int, int)> lambdaInitState = StateInitializer(line, nt, state);
    //- QMDD tree has index 0 at last node. Function works from least significant digit(nt-1) through most significant(0)
    //- lambdas used for converting digit (of input numbber) index to qubit index
    auto c0indice = [nt](int i){ return (nt - 1 - 3 * i);};
    auto a0indice = [nt](int i){ return (nt - 1 - (3 * i + 1));};
    auto b0indice = [nt](int i){ return (nt - 1 - (3 * i + 2));};
    auto c1indice = [nt](int i){ return (nt - 1 - (3 * (i + 1)));};
    //        auto c0indice = [nt](int i){ return (3 * i);};
    //        auto a0indice = [nt](int i){ return ((3 * i + 1));};
    //        auto b0indice = [nt](int i){ return ((3 * i + 2));};
    //        auto c1indice = [nt](int i){ return ((3 * (i + 1)));};
    for(int i = 0; i < n; ++i){
        lambdaInitState(num1[i], a0indice(i));//put first number in register
        lambdaInitState(num2[i], b0indice(i));//put second number in register
    }
    //execute ripple adder:
    auto lambdaCarry = [&state, this, line, nt](int c0, int a0, int b0, int c1){
        gg.ToffoliGenOrApply(line, c1, a0, b0, nt, &state);
        gg.CNotGenOrApply(line, b0, a0, nt, &state);
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
    };
    auto lambdaCarryInv = [&state, this, line, nt](int c0, int a0, int b0, int c1){
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
        gg.CNotGenOrApply(line, b0, a0, nt, &state);
        gg.ToffoliGenOrApply(line, c1, a0, b0, nt, &state);
    };
    auto lambdaSum = [&state, this, line, nt](int c, int a, int b){
        gg.CNotGenOrApply(line, b, a, nt, &state);
        gg.CNotGenOrApply(line, b, c, nt, &state);
    };
    //execute ripple adder:
    for (int i = 0; i < n; ++i){
        lambdaCarry(c0indice(i), a0indice(i), b0indice(i), c1indice(i));
    }
    gg.CNotGenOrApply(line, b0indice(n-1), a0indice(n-1), nt, &state);
    lambdaSum(c0indice(n-1), a0indice(n-1), b0indice(n-1));
    for (int i = n - 2; i >= 0; --i){
        lambdaCarryInv(c0indice(i), a0indice(i), b0indice(i), c1indice(i));
        lambdaSum(c0indice(i), a0indice(i), b0indice(i));
    }
    delete[] line;
    return state;
}
std::function<void (int, int)> RegisterFactory::StateInitializer(short *line, int nt, dd::Edge &state) {
    auto lambdaInitState = [&state, this, line, nt](int val, int index){
        switch (val) {
            case 0:
                break;
            case 1:
                gg.NotGenOrApply(line, index, nt, &state);
                break;
            case 2:
                gg.HadGenOrApply(line, index, nt, &state);
                break;
            case 3:
                gg.NotGenOrApply(line, index, nt, &state);
                gg.HadGenOrApply(line, index, nt, &state);
                break;
            default:
                assert(0);
                break;
        }
    };
    return lambdaInitState;
}

/// RippleAdderHalfClassic, used when one input (a0) is classic. related to figure 2 of paper. state is prepared inside.
/// @param cnum classic number to be added to quantum register.
dd::Edge RegisterFactory::RippleAdderHalfClassicDebug(ulli cnum, vector<int> num){
    assert(n == static_cast<int>(num.size()));
    
    vector<bool> a0clb2 = shor::base2rep(cnum);//classical value in base 2
    int nt = 2 * n + 1;//quantum number (n), carries(n), overflow carry(1).
    short* line = new short[nt];
    gg.lineClear(line, nt);
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    std::function<void (int, int)> lambdaInitState = StateInitializer(line, nt, state);
    //lambdas used for converting digit (of input numbber) index to qubit index
    auto c0indice = [nt](int i){ return nt - 1 - 2 * i;};
    auto a0indice = [a0clb2](int i){return a0clb2.empty() ? 0 : a0clb2[i]; };
    auto b0indice = [nt](int i){ return nt - 1 - (2 * i + 1);};
    auto c1indice = [nt](int i){ return nt - 1 - 2 * (i + 1);};
    for(int i = 0; i < n; ++i){
        lambdaInitState(num[i], b0indice(i));//put quantum number in register
    }
    HelperRippleAdderHalfClassic(a0indice, c0indice, b0indice, c1indice, line, nt, state);
    delete[] line;
    return state;
}
/// Extracted by Xcode, Ripple Adder Half Classic for adder modulo N.
/// @param b0indice b0 indice generated by a lambda
/// @param c0indice c0 indice generated by a lambda
/// @param c1indice c1 indice generated by a lambda
/// @param a0cnumindice a0, classic number indice generated by a lambda
/// @param line 'line' for gate generation
/// @param nt total number of needed variables
/// @param state input state containing input numbers data.
void RegisterFactory::HelperRippleAdderHalfClassic(const std::function<int (int)> & a0cnumindice, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state) {
    
    auto lambdaCarry = [&state, this, line, nt](int c0, int a0, int b0, int c1){
        if(a0){
            gg.CNotGenOrApply(line, c1, b0, nt, &state);
            gg.NotGenOrApply(line, b0, nt, &state);
        }
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
    };
    auto lambdaCarryInv = [&state, this, line, nt](int c0, int a0, int b0, int c1){
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
        if(a0){
            gg.NotGenOrApply(line, b0, nt, &state);
            gg.CNotGenOrApply(line, c1, b0, nt, &state);
        }
    };
    auto lambdaSum = [&state, this, line, nt](int c0, int a0, int b0){
        gg.CNotGenOrApply(line, b0, c0, nt, &state);
        //Problem is here!
        if(a0){
            gg.NotGenOrApply(line, b0, nt, &state);
        }
    };
    //execute ripple adder:
    for (int i = 0; i < n; ++i){
        lambdaCarry(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
    }
    if(a0cnumindice(n-1))
        gg.NotGenOrApply(line, b0indice(n-1), nt, &state);
    
    lambdaSum(c0indice(n-1), a0cnumindice(n-1), b0indice(n-1));
    for (int i = n - 2; i >= 0; --i){
        lambdaCarryInv(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
        lambdaSum(c0indice(i), a0cnumindice(i), b0indice(i));
    }
}
/// Like HelperRippleAdderHalfClassic, int cindex parameter added. qubit cindex negative controls effect of adding classical number to quantum number.
/// @param b0indice b0 indice generated by a lambda
/// @param c0indice c0 indice generated by a lambda
/// @param c1indice c1 indice generated by a lambda
/// @param a0cnumindice a0, classic number indice generated by a lambda
/// @param cindex extra control index
/// @param line 'line' for gate generation
/// @param n num variables to hold input numbers
/// @param nt total number of needed variables
/// @param state input state containing input numbers data.
void RegisterFactory::CRippleAdderHalfClassic(const std::function<int (int)> &a0cnumindice, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice,  int cindex, short *line, int n, int nt, dd::Edge &state){
    auto lambdaCarry = [&state, this, line, nt, cindex](int c0, int a0, int b0, int c1){
        if(a0){
            gg.ToffoliGenOrApply(line, c1, b0, cindex, nt, &state, true, false);
            gg.CNotGenOrApply(line, b0, cindex, nt, &state, false);
        }
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
    };
    auto lambdaCarryInv = [&state, this, line, cindex, nt](int c0, int a0, int b0, int c1){
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
        if(a0){
            gg.CNotGenOrApply(line, b0, cindex, nt, &state, false);
            gg.ToffoliGenOrApply(line, c1, b0, cindex, nt, &state,true ,false);
        }
    };
    auto lambdaSum = [&state, this, line, cindex, nt](int c0, int a0, int b0){
        gg.CNotGenOrApply(line, b0, c0, nt, &state);
        if(a0){
            gg.CNotGenOrApply(line, b0, cindex, nt, &state, false);
        }
    };
    //execute ripple adder:
    for (int i = 0; i < n; ++i){
        lambdaCarry(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
    }
    if(a0cnumindice(n-1))
        gg.CNotGenOrApply(line, b0indice(n-1), cindex, nt, &state, false);
    lambdaSum(c0indice(n-1), a0cnumindice(n-1), b0indice(n-1));
    for (int i = n - 2; i >= 0; --i){
        lambdaCarryInv(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
        lambdaSum(c0indice(i), a0cnumindice(i), b0indice(i));
    }
}
void RegisterFactory::HelperModuloNAdderHalfClassic(const std::function<int (int)> &a0Nindice, const std::function<int (int)> &a0cnumindice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state, const std::function<int ()> &tindice) {
    HelperRippleAdderHalfClassic(a0cnumindice, c0indice, b0indice, c1indice, line, nt, state);
    InvHelperRippleAdderHalfClassic(a0Nindice, c0indice, b0indice, c1indice, line, nt, state);
    gg.NotGenOrApply(line, c1indice(n - 1)/*nt - 2*/, nt, &state);
    gg.CNotGenOrApply(line, tindice(), c1indice(n - 1)/*nt - 2*/, nt, &state);
    gg.NotGenOrApply(line, c1indice(n - 1)/*nt - 2*/, nt, &state);
    CRippleAdderHalfClassic(a0Nindice, c0indice, b0indice, c1indice, tindice(), line, n, nt, state);
    InvHelperRippleAdderHalfClassic(a0cnumindice, c0indice, b0indice, c1indice, line, nt, state);
    gg.CNotGenOrApply(line, tindice(), c1indice(n - 1), nt, &state);
    HelperRippleAdderHalfClassic(a0cnumindice, c0indice, b0indice, c1indice, line, nt, state);
}
/// Controlled-Controlled ModuloNAdder classic number to quantum number.
/// @param mcindex multiplier controller (c in fig 5)
/// @param xindex x index (from quantum number)
/// @param tindex temp memory index
/// @param a0cnumbase2 a(2^i)
/// @param a0Nbase2 N (to be factorized) in base 2
/// @param c0indice current carry
/// @param b0indice current second number in carry element
/// @param c1indice next carry
/// @param line line array for making gates
/// @param nt total number of qubits
/// @param state state to operate on
void RegisterFactory::CCModuloNAdderHalfClassic(int mcindex, int xindex, int tindex, const std::function<int (int)> &a0cnumbase2, const std::function<int (int)> &a0Nbase2, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state){
    CCRippleAdderHalfClassic(mcindex, xindex, tindex, a0cnumbase2, c0indice, b0indice, c1indice, line, nt, state);
    //TODO: take out first subtractor and second adder by introducting CCC adder and multiplier.
    InvHelperRippleAdderHalfClassic(a0Nbase2, c0indice, b0indice, c1indice, line, nt, state);
    gg.NotGenOrApply(line, c1indice(n - 1), nt, &state);
    gg.CNotGenOrApply(line, tindex, c1indice(n - 1), nt, &state);
    gg.NotGenOrApply(line, c1indice(n - 1), nt, &state);
    CRippleAdderHalfClassic(a0Nbase2, c0indice, b0indice, c1indice, tindex, line, n, nt, state);
    InvCCRippleAdderHalfClassic(mcindex, xindex, tindex, a0cnumbase2, c0indice, b0indice, c1indice, line, nt, state);
    gg.CNotGenOrApply(line, tindex, c1indice(n - 1), nt, &state);
    CCRippleAdderHalfClassic(mcindex, xindex, tindex, a0cnumbase2, c0indice, b0indice, c1indice, line, nt, state);
}
/// Controlled-Controlled ModuloNSubtractor classic number from quantum number. Used in exponentiator (fig 6)
/// @param mcindex multiplier controller (c in fig 5)
/// @param xindex x index (from quantum number)
/// @param tindex temp memory index
/// @param a0cnumbase2 a(2^i)
/// @param a0Nbase2 N (to be factorized) in base 2
/// @param c0indice current carry
/// @param b0indice current second number in carry element
/// @param c1indice next carry
/// @param line line array for making gates
/// @param nt total number of qubits
/// @param state state to operate on
void RegisterFactory::InvCCModuloNAdderHalfClassic(int mcindex, int xindex, int tindex, const std::function<int (int)> &a0cnumbase2, const std::function<int (int)> &a0Nbase2, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state){
 InvCCRippleAdderHalfClassic(mcindex, xindex, tindex, a0cnumbase2, c0indice, b0indice, c1indice, line, nt, state);
      gg.CNotGenOrApply(line, tindex, c1indice(n - 1), nt, &state);
      CCRippleAdderHalfClassic(mcindex, xindex, tindex, a0cnumbase2, c0indice, b0indice, c1indice, line, nt, state);
    InvCRippleAdderHalfClassic(a0Nbase2, c0indice, b0indice, c1indice, tindex, line, n, nt, state);
      gg.NotGenOrApply(line, c1indice(n - 1), nt, &state);
      gg.CNotGenOrApply(line, tindex, c1indice(n - 1), nt, &state);
        gg.NotGenOrApply(line, c1indice(n - 1), nt, &state);
      HelperRippleAdderHalfClassic(a0Nbase2, c0indice, b0indice, c1indice, line, nt, state);
        InvCCRippleAdderHalfClassic(mcindex, xindex, tindex, a0cnumbase2, c0indice, b0indice, c1indice, line, nt, state);
}
dd::Edge RegisterFactory::ModuloNAdderHalfClassicDebug(ulli cnum, vector<int> qnum){
    n = static_cast<int>(base2N.size());
    assert(n == qnum.size());//assumed quantum number represented in n base2 digits.
    short* line = new short[n];
    gg.lineClear(line, n);
    vector<bool> a0clb2 = shor::base2rep(cnum, n);//classical value in base 2
    int nt = 2 * n + 2;//quantum number (n), carries(n), overflow carrie(1), temp memory qubit(1).
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    std::function<void (int, int)> lambdaInitState = StateInitializer(line, nt, state);
    //lambdas used for converting digit (of input numbber) index to qubit index
    //a0 line does not exist. replaced by classic register.
    auto c0indice = [nt](int i){ return nt - 2 - 2 * i;};
    auto a0cnumindice = [a0clb2](int i){return a0clb2.empty() ? 0 : a0clb2[i]; };//return cnum in base 2.
    auto a0Nindice = [Nrep = this->base2N](int i){return Nrep.empty() ? 0 : Nrep[i]; };//return cnum in base 2.
    auto b0indice = [nt](int i){ return nt - 2 - (2 * i + 1);};
    auto c1indice = [nt](int i){ return nt - 2 - 2 * (i + 1);};
    auto tindice = [nt](){return nt - 1;};//temp memory of modulon qubit(t in fig 4 of paper)
    for(int i = 0; i < n; ++i){
        lambdaInitState(qnum[i], b0indice(i));//put quantum number in register
    }
    HelperModuloNAdderHalfClassic(a0Nindice, a0cnumindice, b0indice, c0indice, c1indice, line, nt, state, tindice);
    
    delete[] line;
    return state;
}
void RegisterFactory::InvHelperRippleAdderHalfClassic(const std::function<int (int)> &a0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state) {
    auto lambdaCarry = [&state, this, line, nt](int c0, int a0, int b0, int c1){
        if(a0){
            gg.CNotGenOrApply(line, c1, b0, nt, &state);
            gg.NotGenOrApply(line, b0, nt, &state);
        }
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
    };
    auto lambdaCarryInv = [&state, this, line, nt](int c0, int a0, int b0, int c1){
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
        
        if(a0){
            gg.NotGenOrApply(line, b0, nt, &state);
            gg.CNotGenOrApply(line, c1, b0, nt, &state);
        }
    };
    auto lambdaSum = [&state, this, line, nt](int c0, int a0, int b0){
        gg.CNotGenOrApply(line, b0, c0, nt, &state);
        if(a0)
            gg.NotGenOrApply(line, b0, nt, &state);
    };
    for (int i = 0; i < n - 1; ++i){
        lambdaSum(c0indice(i), a0indice(i), b0indice(i));
        lambdaCarry(c0indice(i), a0indice(i), b0indice(i), c1indice(i));
    }
    lambdaSum(c0indice(n - 1), a0indice(n - 1), b0indice(n - 1));
    if(a0indice(n-1))
        gg.NotGenOrApply(line, b0indice(n-1), nt, &state);
    for (int i = n - 1; i >= 0; --i){
        lambdaCarryInv(c0indice(i), a0indice(i), b0indice(i), c1indice(i));
    }
}
/// Controlled Half Classic Multiplier (fig 5)
/// @param a0Nbase2 N (to be factorized) in base 2
/// @param b0indice  quantum register. current second number in carry element
/// @param c0indice current carry
/// @param c1indice next carry
/// @param a0cnum classical number to be multiplied by the quantum number.
/// @param line line array for making gates
/// @param mcindex  multiplier controller (c in fig 5)
/// @param nt total number of qubits
/// @param state state to operate on
/// @param tindex temp memory index for moduloNAdder (FIG2)
/// @param xindice  quantum register. x indice (for input quantum number)
void RegisterFactory::CMultiplierModuloNHalfClassic(const std::function<int (int)> &a0Nbase2, const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, ulli a0cnum, short *line, int mcindex, int nt, dd::Edge &state, const std::function<int ()> &tindex, const std::function<int (int)> &xindice) {
    
    for(int i = 0; i < n; ++i){//fig5 of paper, repeated shift and apply
        //TODO: make next four opertations not do redundant calculations.
        ulli tempa0cnum = a0cnum;
        tempa0cnum = tempa0cnum << i;
        tempa0cnum = tempa0cnum % N;
        vector<bool> a0clb2 = shor::base2rep(tempa0cnum, n);//classical value in base 2, n digits.
        auto a0cnumbase2 = [&a0clb2](int i){return a0clb2.empty() ? 0 : a0clb2[i]; };//return cnum's digits in base 2.
        CCModuloNAdderHalfClassic(mcindex, xindice(i), tindex(), a0cnumbase2, a0Nbase2, c0indice, b0indice ,c1indice, line, nt, state);
    }
    gg.NotGenOrApply(line, mcindex, nt, &state);
    for(int i = 0; i < n; ++i){
        gg.ToffoliGenOrApply(line, b0indice(i), mcindex, xindice(i), nt, &state);
    }
    gg.NotGenOrApply(line, mcindex, nt, &state);
}

/// Inverse Controlled Half Classic Multiplier. Used in Modulo Exponentiator (fig 6)
/// @param a0Nbase2 N (to be factorized) in base 2
/// @param b0indice  quantum register. current second number in carry element
/// @param c0indice current carry
/// @param c1indice next carry
/// @param a0cnum classical number to be multiplied by the quantum number.
/// @param line line array for making gates
/// @param mcindex  multiplier controller (c in fig 5)
/// @param nt total number of qubits
/// @param state state to operate on
/// @param tindex temp memory index for moduloNAdder (FIG2)
/// @param xindice quantum register. x indice (for input quantum number)
void RegisterFactory::InvCMultiplierModuloNHalfClassic(const std::function<int (int)> &a0Nbase2, const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, ulli a0cnum, short *line, int mcindex, int nt, dd::Edge &state, const std::function<int ()> &tindex, const std::function<int (int)> &xindice){
        gg.NotGenOrApply(line, mcindex, nt, &state);
    for(int i = n - 1 ; i >= n; --i){
         gg.ToffoliGenOrApply(line, b0indice(i), mcindex, xindice(i), nt, &state);
     }
    gg.NotGenOrApply(line, mcindex, nt, &state);
    for(int i = n - 1; i >= 0; --i){//fig5 of paper, repeated shift and apply
         //TODO: make next four opertations not do redundant calculations.
        ulli tempa0cnum = a0cnum;
        tempa0cnum = tempa0cnum << i;
        tempa0cnum = tempa0cnum % N;
        vector<bool> a0clb2 = shor::base2rep(tempa0cnum, n);//classical value in base 2, n digits.
        auto a0cnumbase2 = [&a0clb2](int i){return a0clb2.empty() ? 0 : a0clb2[i]; };//return cnum's digits in base 2.
        InvCCModuloNAdderHalfClassic(mcindex, xindice(i), tindex(), a0cnumbase2, a0Nbase2, c0indice, b0indice ,c1indice, line, nt, state);
    }
}
/// Tester of controlled multiplier modulo N. Fig 5 of paper.
/// @param a0cnumbase10 classical number to be multiplied by quantum number.
/// @param qnum quantum number
/// @param mcv multiplier controller qubit value. Switch like in state initializer switch.
dd::Edge RegisterFactory::CMultiplierModuloNDebug(ulli a0cnumbase10, vector<int> qnum, int mcv){
    assert(n == qnum.size());//assert quantum number represented in base2 is n digits.
    int nt = 3 * n + 3;//quantum number x Fig5 (n), quantum register y Fig5 (n), carries(n) Fig2, overflow carry(1) Fig2, temp memory qubit Fig4 t (1), multiplier control qubit c Fig5(1).
    short* line = new short[nt];
    gg.lineClear(line, nt);
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    auto lambdaInitState = StateInitializer(line, nt, state);
    //lambdas used for converting digit (of input number) to qubit index
    //a0 line does not exist. replaced by classic register.
    auto a0Nbase2 = [Nrep = this->base2N](int i){return Nrep.empty() ? 0 : Nrep[i]; };//return N's digits in base 2.
    int mcindex = nt - 1;//multiplier control qubit(c in fig 5 of paper)
    auto xindice = [nt](int i){return nt - 2 - i;};
    auto tindex = [nt, n = this->n](){return nt - n - 2;};
    auto c0indice = [nt, n = this->n](int i){return nt - n - 3 - 2 * i;};
    auto b0indice = [nt, n = this->n](int i){return nt - n - 3 - (2 * i + 1);};
    auto c1indice = [nt, n = this->n](int i){return nt - n - 3 - 2 * (i + 1);};
    
    lambdaInitState(mcv, mcindex);
    for(int i = 0; i < n; ++i){
        lambdaInitState(qnum[i], xindice(i));//put quantum number in register
    }
    
    CMultiplierModuloNHalfClassic(a0Nbase2, b0indice, c0indice, c1indice, a0cnumbase10, line, mcindex, nt, state, tindex, xindice);
    delete[] line;
    return state;
};
/// Helper function extracted by Xcode to be used in CCRippleAdderHalfClassic() and CCRippleSubtractorHalfClassic(). Prepares lambdaCarry, lambdaCarryInv and lambdaSum
/// @param lambdaCarry lambdaCarry
/// @param lambdaCarryInv lambdaCarryInverse
/// @param lambdaSum lambdaSum
/// @param line 'line' for gate generation
/// @param mc multiply control qubit, c in fig 5 of paper
/// @param nt total number of needed variables
/// @param state input state containing input numbers data.
/// @param x individual x qubit index in fig5
void RegisterFactory::HelperCCRippleHalfClassic(std::function<void (int, int, int, int)> &lambdaCarry, std::function<void (int, int, int, int)> &lambdaCarryInv, std::function<void (int, int, int)> &lambdaSum, short *line, int mc, int nt, dd::Edge &state, int x) {
    lambdaCarry = [&state, this, line, mc, x, nt](int c0, int a0, int b0, int c1){
        if(a0){
            map<int, bool> m{{mc, true}, {x, true}, {b0,true}};
            gg.CKNotGenOrApply(line, c1, m, nt, &state);
            gg.ToffoliGenOrApply(line, b0, mc, x, nt, &state);
        }
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
    };
    lambdaCarryInv = [&state, this, line, mc, x, nt](int c0, int a0, int b0, int c1){
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
        if(a0){
            gg.ToffoliGenOrApply(line, b0, mc, x, nt, &state);
            map<int, bool> m{{mc, true}, {x, true}, {b0,true}};
            gg.CKNotGenOrApply(line, c1, m, nt, &state);
        }
    };
    lambdaSum = [&state, this, line, mc, x, nt](int c0, int a0, int b0){
        if(a0){
            gg.ToffoliGenOrApply(line, b0, mc, x, nt, &state);
        }
        gg.CNotGenOrApply(line, b0, c0, nt, &state);
    };
}

/// Used in controlled multiplication modulo N. Fig 5 of paper
/// @param mc multiply control qubit, c in fig 5 of paper
/// @param x individual x qubit index in fig5
/// @param t temp memory in qubit in fig 4
/// @param a0cnumindice a0 (current first number, classic). a2^i in fig 5.
/// @param c0indice c0 (current carry) produced by lambda
/// @param b0indice b0 line (current second number) produced by lambda
/// @param c1indice c1 (next carry) produced by lambda
/// @param line 'line' for gate generation
/// @param nt total number of needed variables
/// @param state input state containing input numbers data.
void RegisterFactory::CCRippleAdderHalfClassic(int mc, int x, int t, const std::function<int (int)> &a0cnumindice, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state){
    std::function<void (int, int, int, int)> lambdaCarry;
    std::function<void (int, int, int, int)> lambdaCarryInv;
    std::function<void (int, int, int)> lambdaSum;
    HelperCCRippleHalfClassic(lambdaCarry, lambdaCarryInv, lambdaSum, line, mc, nt, state, x);
    //execute ripple adder:
    for (int i = 0; i < n; ++i){
        lambdaCarry(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
    }
    if(a0cnumindice(n-1))
        gg.ToffoliGenOrApply(line, b0indice(n-1), mc, x, nt, &state);
    lambdaSum(c0indice(n-1), a0cnumindice(n-1), b0indice(n-1));
    for (int i = n - 2; i >= 0; --i){
        lambdaCarryInv(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
        lambdaSum(c0indice(i), a0cnumindice(i), b0indice(i));
    }
}
/// Used in controlled multiplication modulo N. Fig 5 of paper
/// @param mc multiply control qubit, c in fig 5 of paper
/// @param x x qubits in fig5
/// @param t temp memory in qubit in fig 4
/// @param a0cnumindice a0 (current first number, classic). a2^i in fig 5.
/// @param c0indice c0 (current carry) produced by lambda
/// @param b0indice b0 line (current second number) produced by lambda
/// @param c1indice c1 (next carry) produced by lambda
/// @param line 'line' for gate generation
/// @param nt total number of needed variables
/// @param state input state containing input numbers data.
void RegisterFactory::InvCCRippleAdderHalfClassic(int mc, int xindex, int t, const std::function<int (int)> &a0cnumindice, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state){
    std::function<void (int, int, int, int)> lambdaCarry;
    std::function<void (int, int, int, int)> lambdaCarryInv;
    std::function<void (int, int, int)> lambdaSum;
    HelperCCRippleHalfClassic(lambdaCarry, lambdaCarryInv, lambdaSum, line, mc, nt, state, xindex);
    for (int i = 0; i < n - 1; ++i){
        lambdaSum(c0indice(i), a0cnumindice(i), b0indice(i));
        lambdaCarry(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
    }
    lambdaSum(c0indice(n-1), a0cnumindice(n-1), b0indice(n-1));
    if(a0cnumindice(n-1))
        gg.ToffoliGenOrApply(line, b0indice(n-1), mc, xindex, nt, &state);
    for (int i = n - 1; i >= 0; --i){
        lambdaCarryInv(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
    }
}
/// ExponentiatorModuloNDebugger(tester), Fig6 of the paper.
/// @param cnum classic number to be base for exponentiator. 'a' in fig6 of paper.
/// @param qnum quantum number x in fig6, which will be the exponent. coded with arbitrary qubit size in this function.
dd::Edge RegisterFactory::ExponentiatorModuloNDebug(ulli cnum, vector<int> qnum){
    int m = static_cast<int>(qnum.size());
    //TODO: Map of {i \in 0 \to 2m} to a^2^i and a^2^-i
    int nt = 3 * n + 2 + m;//quantum number x Fig5 (n), quantum register y Fig5 (n), carries(n) Fig2, overflow carry(1) Fig2, temp memory qubit Fig4 t (1), quantum exponent(m) which are controls for multipliers.
    short* line = new short[nt];
    gg.lineClear(line, nt);
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    vector<bool> a0clb2 = shor::base2rep(1, n);//classical value in base 2, n digits.
    auto lambdaInitState = StateInitializer(line, nt, state);
    //a0 line does not exist. replaced by classic register.
    auto a0Nbase2 = [Nrep = this->base2N](int i){return Nrep.empty() ? 0 : Nrep[i]; };//return N's digits in base 2.
    //lambdas used for converting digit (of input number) to qubit index.
    auto xindice = [nt](int i){return nt - i - 1;};//x qubits in fig 6. quantum register
    auto xpindice = [nt, m](int i){return nt - m - 1 - i;};//= x qubits in fig 5
    auto c0indice = [nt, m, n = this->n](int i){return nt - m - n - 1 - 2 * i;};
    auto b0indice = [nt, m, n = this->n](int i){return nt - m - n - 1 - (2 * i + 1);}; //quantum register
    auto c1indice = [nt, m, n = this->n](int i){return nt - m - n - 1 - 2 * (i + 1);};
    auto tindex = [nt, m, n = this->n](){return nt - m - 3 * n - 2;};
    for(int i = 0; i < n; ++i){
           lambdaInitState(qnum[i], xindice(i));//put quantum number in register
       }
    gg.NotGenOrApply(line, xpindice(n-1), nt);//put 1 in result register.
    ulli clop = 1;//classical operand of multipliers
    for(int i = 0; i < m; ++i){//repeatedly multiply a^2^i and a^2^-i, based on control qubit (x_i) as in fig 6 of paper.
        CMultiplierModuloNHalfClassic(a0Nbase2, b0indice, c0indice, c1indice, clop, line, xindice(i), nt, state, tindex, xpindice);
        ulli invclop = shor::modInverse(clop, N);
        InvCMultiplierModuloNHalfClassic(a0Nbase2, xpindice, c0indice, c1indice, invclop, line, xindice(i), nt, state, tindex, b0indice);
        clop = clop * clop;
    }
    delete[] line;
    return state;
}
RegisterFactory::~RegisterFactory(){
    delete[] line;
}
/// inverse of CRippleAdderHalfClassic used in exponentiator, fig6
/// @param b0indice b0 indice generated by a lambda
/// @param c0indice c0 indice generated by a lambda
/// @param c1indice c1 indice generated by a lambda
/// @param a0cnumindice a0, classic number indice generated by a lambda
/// @param cindex extra control index
/// @param line 'line' for gate generation
/// @param n num variables to hold input numbers
/// @param nt total number of needed variables
/// @param state input state containing input numbers data.
void RegisterFactory::InvCRippleAdderHalfClassic(const std::function<int (int)> &a0cnumindice, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice,  int cindex, short *line, int n, int nt, dd::Edge &state) {
    auto lambdaCarry = [&state, this, line, nt, cindex](int c0, int a0, int b0, int c1){
        if(a0){
            gg.ToffoliGenOrApply(line, c1, b0, cindex, nt, &state, true, false);
            gg.CNotGenOrApply(line, b0, cindex, nt, &state, false);
        }
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
    };
    auto lambdaCarryInv = [&state, this, line, cindex, nt](int c0, int a0, int b0, int c1){
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
        if(a0){
            gg.CNotGenOrApply(line, b0, cindex, nt, &state, false);
            gg.ToffoliGenOrApply(line, c1, b0, cindex, nt, &state,true ,false);
        }
    };
    auto lambdaSum = [&state, this, line, cindex, nt](int c0, int a0, int b0){
         gg.CNotGenOrApply(line, b0, c0, nt, &state);
         if(a0){
             gg.CNotGenOrApply(line, b0, cindex, nt, &state, false);
         }
     };
        for (int i = n - 2; i >= 0; --i){
                    lambdaSum(c0indice(i), a0cnumindice(i), b0indice(i));
        lambdaCarry(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
    };
    lambdaSum(c0indice(n-1), a0cnumindice(n-1), b0indice(n-1));
    if(a0cnumindice(n-1))
    gg.CNotGenOrApply(line, b0indice(n-1), cindex, nt, &state, false);
      for (int i = 0; i < n; ++i){
          lambdaCarryInv(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
      }
}

