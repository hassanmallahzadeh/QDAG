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
    std::function<void (int, int)> lambdaInitState = ExtractedStateInitializer(line, nt, state);
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
std::function<void (int, int)> RegisterFactory::ExtractedStateInitializer(short *line, int nt, dd::Edge &state) {
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
    std::function<void (int, int)> lambdaInitState = ExtractedStateInitializer(line, nt, state);
    //lambdas used for converting digit (of input numbber) index to qubit index
    auto c0indice = [nt](int i){ return nt - 1 - 2 * i;};
    auto a0indice = [a0clb2](int i){return a0clb2.empty() ? 0 : a0clb2[i]; };
    auto b0indice = [nt](int i){ return nt - 1 - (2 * i + 1);};
    auto c1indice = [nt](int i){ return nt - 1 - 2 * (i + 1);};
    
    for(int i = 0; i < n; ++i){
        lambdaInitState(num[i], b0indice(i));//put quantum number in register
    }
    HelperRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0indice, line, n, nt, state);
    delete[] line;
    return state;
}
/// Extracted by Xcode, Ripple Adder Half Classic for adder modulo N.
/// @param b0indice b0 indice generated by a lambda
/// @param c0indice c0 indice generated by a lambda
/// @param c1indice c1 indice generated by a lambda
/// @param a0cnumindice a0, classic number indice generated by a lambda
/// @param line 'line' for gate generation
/// @param n num variables to hold input numbers
/// @param nt total number of needed variables
/// @param state input state containing input numbers data.
void RegisterFactory::HelperRippleAdderHalfClassic(const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, const std::function<int (int)> & a0cnumindice, short *line, int n, int nt, dd::Edge &state) {
    
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
///// Like ExtractedRippleAdderHalfClassicV1, int t parameter added. qubit t negative controls effect of adding classical number to quantum number.
void RegisterFactory::CRippleAdderHalfClassic(const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, const std::function<int (int)> &a0indice, int t, short *line, int n, int nt, dd::Edge &state){
    auto lambdaCarry = [&state, this, line, nt, t](int c0, int a0, int b0, int c1){
        if(a0){
            gg.ToffoliGenOrApply(line, c1, b0, t, nt, &state, true, false);
            gg.CNotGenOrApply(line, b0, t, nt, &state, false);
        }
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
    };
    auto lambdaCarryInv = [&state, this, line, t, nt](int c0, int a0, int b0, int c1){
        gg.ToffoliGenOrApply(line, c1, c0, b0, nt, &state);
        if(a0){
            gg.CNotGenOrApply(line, b0, t, nt, &state, false);
            gg.ToffoliGenOrApply(line, c1, b0, t, nt, &state,true ,false);
        }
    };
    auto lambdaSum = [&state, this, line, t, nt](int c0, int a0, int b0){
        gg.CNotGenOrApply(line, b0, c0, nt, &state);
        if(a0){
            gg.CNotGenOrApply(line, b0, t, nt, &state, false);
        }
    };
    //execute ripple adder:
    for (int i = 0; i < n; ++i){
        lambdaCarry(c0indice(i), a0indice(i), b0indice(i), c1indice(i));
    }
    if(a0indice(n-1))
        gg.CNotGenOrApply(line, b0indice(n-1), t, nt, &state, false);
    lambdaSum(c0indice(n-1), a0indice(n-1), b0indice(n-1));
    for (int i = n - 2; i >= 0; --i){
        lambdaCarryInv(c0indice(i), a0indice(i), b0indice(i), c1indice(i));
        lambdaSum(c0indice(i), a0indice(i), b0indice(i));
    }
}
void RegisterFactory::HelperModuloNAdderHalfClassic(const std::function<int (int)> &a0Nindice, const std::function<int (int)> &a0cnumindice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state, const std::function<int ()> &tindice) {
    HelperRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0cnumindice, line, n, nt, state);
    HelperRippleSubtractorHalfClassic(b0indice, c0indice, c1indice, a0Nindice, line, nt, state);
    gg.NotGenOrApply(line, c1indice(n - 1)/*nt - 2*/, nt, &state);
    gg.CNotGenOrApply(line, tindice(), c1indice(n - 1)/*nt - 2*/, nt, &state);
    gg.NotGenOrApply(line, c1indice(n - 1)/*nt - 2*/, nt, &state);
    CRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0Nindice, tindice(), line, n, nt, state);
    HelperRippleSubtractorHalfClassic(b0indice, c0indice, c1indice, a0cnumindice, line, nt, state);
    gg.CNotGenOrApply(line, tindice(), c1indice(n - 1), nt, &state);
    HelperRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0cnumindice, line, n, nt, state);
}
void RegisterFactory::CCModuloNAdderHalfClassic(int mcindex, int xindex, int tindex, const std::function<int (int)> &a0cnumindice, const std::function<int (int)> &a0Nindice, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state){
    CCRippleAdderHalfClassic(mcindex, xindex, tindex, a0cnumindice, c0indice, b0indice, c1indice, line, nt, state);
    HelperRippleSubtractorHalfClassic(a0Nindice, c0indice, b0indice, c1indice, line, nt, state);
    gg.NotGenOrApply(line, c1indice(n - 1), nt, &state);
    gg.CNotGenOrApply(line, tindex, c1indice(n - 1), nt, &state);
    gg.NotGenOrApply(line, c1indice(n - 1), nt, &state);
    CRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0Nindice, tindex, line, n, nt, state);
    CCRippleSubtractorHalfClassic(mcindex, xindex, tindex, a0cnumindice, c0indice, b0indice, c1indice, line, nt, state);
    gg.CNotGenOrApply(line, tindex, c1indice(n - 1), nt, &state);
    CCRippleAdderHalfClassic(mcindex, xindex, tindex, a0cnumindice, c0indice, b0indice, c1indice, line, nt, state);
}

dd::Edge RegisterFactory::ModuloNAdderHalfClassicDebug(ulli cnum, vector<int> qnum){
    n = static_cast<int>(base2N.size());
    assert(n == qnum.size());//assumed quantum number represented in n base2 digits.
    short* line = new short[n];
    gg.lineClear(line, n);
    vector<bool> a0clb2 = shor::base2rep(cnum, n);//classical value in base 2
    int nt = 2 * n + 2;//quantum number (n), carries(n), overflow carrie(1), temp memory qubit(1).
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    std::function<void (int, int)> lambdaInitState = ExtractedStateInitializer(line, nt, state);
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
void RegisterFactory::HelperRippleSubtractorHalfClassic(const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, const std::function<int (int)> &a0indice, short *line, int nt, dd::Edge &state) {
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
dd::Edge RegisterFactory::CMultiplierModuloNClassicDebug(ulli cnum, vector<int> qnum){
    n = static_cast<int>(base2N.size());
    assert(n == qnum.size());//assumed quantum number represented in n base2 digits.
    short* line = new short[n];
    gg.lineClear(line, n);
    vector<bool> a0clb2 = shor::base2rep(cnum, n);//classical value in base 2
    int nt = 2 * n + 3;//quantum number (n), carries(n), overflow carrie(1), temp memory qubit(1), multiplier control qubit(1).
    
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    auto lambdaInitState = ExtractedStateInitializer(line, nt, state);
    //lambdas used for converting digit (of input numbber) index to qubit index
    //a0 line does not exist. replaced by classic register, used in moduloN circuit.
     auto a0Nindice = [Nrep = this->base2N](int i){return Nrep.empty() ? 0 : Nrep[i]; };//return cnum in base 2.
    auto mcindex = [nt](){return nt - 1;};//multiplier control qubit(c in fig 5 of paper)
    auto xindice = [nt](int i){ return nt - 2 - i;};
    auto tindex = [nt, n = this->n](){return nt - n - 1;};
    auto c0indice = [nt, n = this->n](int i){ return nt - n - 2 - 2 * i;};
    auto a0cnumindice = [a0clb2](int i){return a0clb2.empty() ? 0 : a0clb2[i]; };//return cnum in base 2.
    auto b0indice = [nt, n = this->n](int i){ return nt - n - 2 - (2 * i + 1);};
    auto c1indice = [nt, n = this->n](int i){ return nt - n - 2 - 2 * (i + 1);};
    
    for(int i = 0; i < n; ++i){
        lambdaInitState(qnum[i], xindice(i));//put quantum number in register
    }
    for(int i = 0; i < n; ++i){//fig5 of paper, repeated shift and apply
        CCModuloNAdderHalfClassic(mcindex(), xindice(i), tindex(),a0cnumindice, a0Nindice, c0indice ,b0indice ,c1indice ,line, nt, state);
        cnum <<= 1;
        a0clb2 = shor::base2rep(cnum, n);//TODO: look how to do this more efficiently. We don't have to calculate the bit pattern again after a shift.
    }
    gg.NotGenOrApply(line, mcindex(), nt);
    for(int i = 0; i < n; ++i){
        gg.ToffoliGenOrApply(line, b0indice(i), mcindex(), xindice(i), nt);
    }
    gg.NotGenOrApply(line, mcindex(), nt);
    return state;
};
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
/// @param x x qubits in fig5
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
void RegisterFactory::CCRippleSubtractorHalfClassic(int mc, int x, int t, const std::function<int (int)> &a0cnumindice, const std::function<int (int)> &c0indice, const std::function<int (int)> &b0indice, const std::function<int (int)> &c1indice, short *line, int nt, dd::Edge &state){
    std::function<void (int, int, int, int)> lambdaCarry;
       std::function<void (int, int, int, int)> lambdaCarryInv;
       std::function<void (int, int, int)> lambdaSum;
       HelperCCRippleHalfClassic(lambdaCarry, lambdaCarryInv, lambdaSum, line, mc, nt, state, x);
    for (int i = 0; i < n - 1; ++i){
                  lambdaSum(c0indice(i), a0cnumindice(i), b0indice(i));
             lambdaCarry(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
         }
    lambdaSum(c0indice(n-1), a0cnumindice(n-1), b0indice(n-1));
      if(a0cnumindice(n-1))
            gg.ToffoliGenOrApply(line, b0indice(n-1), mc, x, nt, &state);
    for (int i = n - 1; i >= 0; --i){
           lambdaCarry(c0indice(i), a0cnumindice(i), b0indice(i), c1indice(i));
       }
}
RegisterFactory::~RegisterFactory(){
    delete[] line;
}
