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
RegisterFactory::RegisterFactory(lli N, lli a, dd::Package *dd) : gg(dd) {
    this->dd = dd;
    if(N > 0){//0 blocks for debugging of ripple carry
        this->N = N;
        this->a = a;
        base2N = shor::base2rep(N);
        n = base2N.size();
        //    DetermineRegSizes();
        //    AllocateRegVectors();
        //  MakeInitialState();
    }
}
/// CExponentiation the input register, put in output register.
dd::Edge RegisterFactory::CExponentiation(){
    dd::Edge e;
    lli cf = a;//constant factor
    lli cfa = 0;//classical factor argument in multiplication
    for(int i = 0; i < n_in; ++i){
        if(i == 0){
            cfa = 1;
        }
        else{
            cf = (cf * cf) % N;
            cfa = cf;
        }
        CMulModN(v_in[i], cfa);
        gg.swapRegistersOnState(n_tot, v_out, v_mod_temp, state);
        CMulModNInv(v_in[i], cfa);
    }
    return e;
}

/// Make initial state by assuming order: input register, then auxiliary register, then output register.
//void RegisterFactory::MakeInitialState(){
//    state = StateGenerator(dd).dd_BaseState(n_tot, 1);//put output reg in
//    line = new short[n_tot]{};// set 'line' for dd->makeGateDD
//    for (int i = 0; i < n_tot; ++i){
//        line[i] = -1;
//    }
//     lli tempN = N;//put N in
//        int i = 0;
//        while(tempN > 0){
//            if(tempN & 1){
//                gg.lineSet(line, n_addmod - i - 1);
//                state = dd->multiply(dd->makeGateDD(Xmat, n_tot, line) , state);
//                gg.lineReset(line, n_addmod - i - 1);
//            }
//            ++i;
//            tempN >>= 1;
//        }
//    assert(i < n_addmod);
//    for(int i = 0; i < n_in; ++i){//prepare input register.
//        gg.lineSet(line, n_out + n_addmod + i);
//        state = dd->multiply(dd->makeGateDD(Hmat, n_tot, line) , state);
//        gg.lineReset(line, n_out + n_addmod + i);
//    }
//};
/// i is the smallest power of 2 needed to store N. Fig 2 of paper.
void RegisterFactory::DetermineRegSizes(){
    lli p = 1;
    int i = 0;//bit count
    do{
        p*=2;
        ++i;
    } while(p < N);
    
    n_out = i;
    n_mul = n_out;
    n_in = 2 * n_out;
    n_car = n_out + 1;
    n_addmod = n_out + 1;//1 to accomodate for single overflow bit in Adder modulo N
    n_mod = n_out;// This register is shared between, modular exponetiator, modular multiplier and modular adder.
    n_tot = n_in + n_car + n_out + n_out + n_addmod + n_mod;
    
}
/// this is to be modified for more efficiency. I am coding the most straightforward implementation
void RegisterFactory::AllocateRegVectors(){//
    for(int i = 0; i < n_addmod; ++i){//first temp register holding N in modulo adder.
        v_addmod_temp.push_back(i);
    }
    for(int i = 0; i < n_mul; ++i){//then temp register holding temp value for multiplication
        v_mul_temp.push_back(n_addmod + i);
    }
    for(int i = 0; i < n_in; ++i){//then input register.
        v_in.push_back(n_mul + n_addmod + i);
    }
    for (int i = 0; i < n_car; ++i){//then aux carry register
        v_car.push_back(n_mul + n_addmod + n_in + i);
    }
    for (int i = 0; i < n_mod; ++i){//then aux temp exponentiation register
        v_mod_temp.push_back(n_mul + n_addmod + n_in + n_car + i);
    }
    for(int i = 0; i < n_out; ++i){//then output register.
        v_out.push_back(n_mul + n_addmod + n_in + n_car + n_mod + i);
    }
}
/// Carry computer gate
/// @param c0 barry from former block
/// @param a0 first non-carry input
/// @param b0 second non-carry input
/// @param c1 current cary
void RegisterFactory::Carry(int c0, int a0, int b0, int c1){
    dd::Edge gate = gg.ToffoliGenOrApply(line, c1, a0, b0, n_tot);
    state = dd->multiply(gate, state);
    gate = gg.CNotGenOrApply(line, b0, a0, n_tot);
    state = dd->multiply(gate, state);
    gate = gg.ToffoliGenOrApply(line, c1, c0, b0, n_tot);
    state = dd->multiply(gate, state);
}
/// Carry Inverse computer gate
/// @param c0 barry from former block
/// @param a0 first non-carry input
/// @param b0 second non-carry input
/// @param c1 current cary
void RegisterFactory::CarryInv(int c0, int a0, int b0, int c1){
    dd::Edge gate = gg.ToffoliGenOrApply(line, c1, c0, b0, n_tot);
    state = dd->multiply(gate, state);
    gate = gg.CNotGenOrApply(line, b0, a0, n_tot);
    state = dd->multiply(gate, state);
    gate = gg.ToffoliGenOrApply(line, c1, a0, b0, n_tot);
    state = dd->multiply(gate, state);
}
/// Sum gate
/// @param c barry from former block
/// @param a first non-carry input
/// @param b second non-carry input (and 'main' output)
void RegisterFactory::Sum(int c, int a, int b){
    dd::Edge gate = gg.CNotGenOrApply(line, b, a, n_tot);
    state = dd->multiply(gate, state);
    gate = gg.CNotGenOrApply(line, b, c, n_tot);
    state = dd->multiply(gate, state);
}
/// Ripple carry adder. Most straight forward implementation based on fig 2 of paper. last line is c_n+1 instead of b_n+1
/// @param v1 first group of qubits (register) for ripple add
/// @param v2 second group of qubits (register) for ripple add
void RegisterFactory::RippleAdder(vector<int> v1, vector<int> v2){
    for(int i = 0; i < n_out; ++i){
        Carry(v_car[i], v1[i], v2[i], v_car[i + 1]);
    }
    dd::Edge gate = gg.CNotGenOrApply(line, v1[v1.size()-1], v2[v2.size()-1], n_tot);
    state = dd->multiply(gate, state);
    Sum(v_car[v_car.size() - 2], v1[v1.size()-1], v2[v2.size()-1]);
    for(int i = n_out - 2 ; i >= 0 ; --i){
        CarryInv(v_car[i], v1[i], v2[i], v_car[i + 1]);
        Sum(v_car[i], v1[i], v2[i]);
    }
}
/// Inverse ripple carry adder. Most straight forward implementation based on fig 2 of paper. last line is c_n+1 instead of b_n+1
/// @param v1 first group of qubits (register) for ripple add
/// @param v2 second group of qubits (register) for ripple add
void RegisterFactory::RippleAdderInv(vector<int> v1, vector<int> v2){
    for(int i = 0 ; i < n_out - 1; ++i){
        Sum(v_car[i], v1[i], v2[i]);
        Carry(v_car[i], v1[i], v2[i], v_car[i + 1]);
    }
    Sum(v_car[v_car.size() - 2], v1[v1.size()-1], v2[v2.size()-1]);
    dd::Edge gate = gg.CNotGenOrApply(line, v1[v1.size()-1], v2[v2.size()-1], n_tot);
    state = dd->multiply(gate, state);
    for(int i = n_out - 1; i >= 0; --i){
        CarryInv(v_car[i], v1[i], v2[i], v_car[i + 1]);
    }
}
/// ModuloNAdder, figure 4 of the paper.
/// @param v1 first group of qubits for inverse add module
/// @param v2  second group of qubits for inverse add module
void RegisterFactory::ModuloNAdder(vector<int> v1, vector<int> v2){
    RippleAdder(v1, v2);
    int temp = v_addmod_temp[v_addmod_temp.size()-1];//last element is temp qubit of ModuloNAdder, does not participate in swap.
    v_addmod_temp.pop_back();
    gg.swapRegistersOnState(n_tot, v1, v_addmod_temp, state);
    v_addmod_temp.push_back(temp);
    RippleAdderInv(v1, v2);
    dd::Edge gate0 = gg.NotGenOrApply(line, v2[v2.size() - 1], n_tot);
    dd::Edge gate1 = gg.CNotGenOrApply(line, v_addmod_temp[n_addmod - 1], v2[v2.size() - 1], n_tot);
    state = dd->multiply(gate0, state);
    state = dd->multiply(gate1, state);
    state = dd->multiply(gate0, state);
    TempNResetter();
    RippleAdder(v1,v2);
    TempNResetter();
    temp = v_addmod_temp[v_addmod_temp.size()-1];//last element is temp qubit of ModuloNAdder, does not participate in swap.
    v_addmod_temp.pop_back();
    gg.swapRegistersOnState(n_tot, v1, v_addmod_temp, state);
    v_addmod_temp.push_back(temp);
    RippleAdderInv(v1, v2);
    state = dd->multiply(gate1, state);
    RippleAdder(v1,v2);
}
/// Inverse ModuloNAdder, figure 4 of the paper.
/// @param v1 first group of qubits for inverse add module
/// @param v2  second group of qubits for inverse add module
void RegisterFactory::ModuloNAdderInv(vector<int> v1, vector<int> v2){
    RippleAdder(v1,v2);
    dd::Edge gate1 = gg.CNotGenOrApply(line, v_addmod_temp[n_addmod - 1], v2[v2.size() - 1], n_tot);
    state = dd->multiply(gate1, state);
    RippleAdderInv(v1, v2);
    int temp = v_addmod_temp[v_addmod_temp.size()-1];//last element of vn is temp qubit of ModuloNAdder
    v_addmod_temp.pop_back();
    gg.swapRegistersOnState(n_tot, v1, v_addmod_temp, state);
    v_addmod_temp.push_back(temp);
    TempNResetter();
    RippleAdder(v1,v2);
    TempNResetter();
    dd::Edge gate0 = gg.NotGenOrApply(line, v2[v2.size() - 1], n_tot);
    state = dd->multiply(gate0, state);
    state = dd->multiply(gate1, state);
    state = dd->multiply(gate0, state);
    RippleAdderInv(v1, v2);
    temp = v_addmod_temp[v_addmod_temp.size()-1];//last element of vn is temp qubit of ModuloNAdder
    v_addmod_temp.pop_back();
    gg.swapRegistersOnState(n_tot, v1, v_addmod_temp, state);
    v_addmod_temp.push_back(temp);
}
/// Function used to set temporaty register of ModuloNAdder to zero (down from N)
void RegisterFactory::TempNResetter(){
    lli tempN = N;
    int i = 0;
    while(tempN > 0){
        if(tempN & 1){
            dd::Edge gate = gg.CNotGenOrApply(line, i, v_addmod_temp[n_addmod-1], n_tot);
            state = dd->multiply(gate, state);
        }
        ++i;
        tempN >>= 1;
    }
    assert(i <= n_out);
}
/// Map a classical value to a register, conditional on controls c1 and c2.
/// @param v target vector for mapping
/// @param c1 control index 1
/// @param cl classical value to be mapped.
/// @param c2 control index 2.
void RegisterFactory::extractedMapToTemp(vector<int> v, int c1, int c2, lli cl) {
    int j = 0;
    while(cl > 0){
        bool apply = cl & 1;
        cl >>= 1;
        if(apply){
            dd::Edge gate =  gg.ToffoliGenOrApply(line, v_mul_temp[j], c1, c2, n_tot);
            state = dd->multiply(gate, state);
        }
        ++j;
        assert(j <= n_out/*or n_mod - 1*/);
    }
}

/// Controlled Multiplication Mod N.
/// @param cq control qubit index. same as power value (i) in a^2^i
/// @param cfa  classical factor argument in multiplication
void RegisterFactory::CMulModN(int cq, lli cfa){
    for(int i = 0; i < n_out; ++i){
        extractedMapToTemp(v_mul_temp, cq, v_out[i], cfa);
        ModuloNAdder(v_mul_temp, v_mod_temp);
        extractedMapToTemp(v_mul_temp, cq, v_out[i], cfa);
        cfa = (cfa * 2) % N;//load the temp register with 'repeated shifted' value
    }
    dd::Edge notGate = gg.NotGenOrApply(line, cq, n_tot);
    state = dd->multiply(notGate, state);
    for(int i = 0; i < n_out; ++i){//copy state from 'x' to 'y' figure 5 of the paper.
        dd::Edge gate =  gg.ToffoliGenOrApply(line, v_mod_temp[i], cq, v_out[i], n_tot);
        state = dd->multiply(gate, state);
    }
    state = dd->multiply(notGate, state);
}
/// Inverse of Controlled Multiplication Mod N.
/// @param c control qubit index
/// @param cfa classical factor argument in multiplication
void RegisterFactory::CMulModNInv(int c, lli cfa){
    dd::Edge notGate = gg.NotGenOrApply(line, c, n_tot);
    state = dd->multiply(notGate, state);
    for(int i = 0; i < n_out; ++i){//copy state from 'x' to 'y' figure 5 of the paper.
        dd::Edge gate =  gg.ToffoliGenOrApply(line, v_mod_temp[i], c, v_out[i], n_tot);//I am assuming Toffolis with different target commute.
        state = dd->multiply(gate, state);
    }
    state = dd->multiply(notGate, state);
    lli cfainv = shor::modInverse(cfa, N);
    vector<lli> v;
    v.push_back(cfainv);
    for (int i = 1; i < n_out; ++i){
        cfainv = (cfainv * 2) % N;
        v.push_back(cfainv);
    }
    for(int i = n_out - 1; i >= 0; --i){
        //load the temp register with (2^i)a
        extractedMapToTemp(v_mul_temp, c, v_out[i], v[i]);
        ModuloNAdderInv(v_mul_temp, v_mod_temp);
        extractedMapToTemp(v_mul_temp, c, v_out[i], v[i]);
    }
}
/// Used to see if ripple adder works as expected. Based on figure 2 of the paper.
///@param num1 first number in addition. size is base 2 representation. for each element, 0 means variable = 0, 1 means varibale = 1, 2 means (|0> + |1>)√2, 3 means |0> - |1>)√2
///@param num2 second number in addition, description like num1
dd::Edge RegisterFactory::RippleAdderDebug(vector<int> num1, vector<int> num2){
    assert(num1.size() == num2.size());
    int n = static_cast<int> (num1.size());
    
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

/// full quantum ripple adder.
/// @param state state for RippleAdder circuit to be executed on.
/// @param n number of digits of numbers to be added (base 2)
dd::Edge RegisterFactory::RippleAdderGeneral(dd::Edge state, int n){
    short* line = new short[n];
    gg.lineClear(line, n);
    int nt = n * 3 + 1;//first number, second number, carry and overflow carry.
    assert(log2(shor::bddNumVar(state, true)) + 1 == nt);
    //lambdas used for converting digit (of input numbber) index to qubit index
    auto c0indice = [](int i){ return 3 * i;};
    auto a0indice = [](int i){ return 3 * i + 1;};
    auto b0indice = [](int i){ return 3 * i + 2;};
    auto c1indice = [](int i){ return 3 * (i + 1);};
    
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
    return gg.permuteOperatorOnState(nt, state);
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
dd::Edge RegisterFactory::RippleAdderHalfClassicDebug(lli cnum, vector<int> num){
    n = static_cast<int>(num.size());
    
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
    ExtractedRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0indice, line, n, nt, state);
    delete[] line;
    return state;
}
/// Extracted by Xcode, Ripple Adder Half Classic for adder modulo N.
/// @param b0indice b0 indice generated by a lambda
/// @param c0indice c0 indice generated by a lambda
/// @param c1indice c1 indice generated by a lambda
/// @param a0indice a0 indice generated by a lambda
/// @param line 'line' for gate generation
/// @param n num variables to hold input numbers
/// @param nt total number of needed variables
/// @param state input state containing input numbers data.
void RegisterFactory::ExtractedRippleAdderHalfClassic(const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, const std::function<int (int)> & a0indice, short *line, int n, int nt, dd::Edge &state) {
    
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
        lambdaCarry(c0indice(i), a0indice(i), b0indice(i), c1indice(i));
    }
    if(a0indice(n-1))
        gg.NotGenOrApply(line, b0indice(n-1), nt, &state);
    
    lambdaSum(c0indice(n-1), a0indice(n-1), b0indice(n-1));
    for (int i = n - 2; i >= 0; --i){
        lambdaCarryInv(c0indice(i), a0indice(i), b0indice(i), c1indice(i));
        lambdaSum(c0indice(i), a0indice(i), b0indice(i));
    }
}
///// Like ExtractedRippleAdderHalfClassicV1, int t parameter added. qubit t negative controls effect of adding classical number to quantum number.
void RegisterFactory::ExtractedControlledRippleAdderHalfClassic(const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, const std::function<int (int)> &a0indice, int t, short *line, int n, int nt, dd::Edge &state){
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
/// Ripple adder half classic. State is prepared in caller.
/// @param cnum number to be added
/// @param state quantum state that number is to be added to
dd::Edge RegisterFactory::RippleAdderHalfClassicGeneral(lli cnum, dd::Edge state){
    int nt = n * 2 + 1;//quantum register, carries and overflow qubit
    vector<bool> a0clb2 = shor::base2rep(cnum, n);//classical value in base 2
    assert(shor::bddNumVar(state, true) == nt);
    
    short* line = new short[n];
    gg.lineClear(line, n);
    
    //lambdas used for converting digit (of input numbber) index to qubit index
    // QMDD tree has index 0 at last node. Function works from least significant digit(nt-1) through most significant(0)
    auto c0indice = [nt](int i){ return nt - 1 - 2 * i;};
    auto a0incide = [a0clb2](int i){return  a0clb2.empty() ? 0 : a0clb2[i];};
    auto b0indice = [nt](int i){ return nt - 1 - (2 * i + 1);};
    auto c1indice = [nt](int i){ return nt - 1 - 2 * (i + 1);};
    ExtractedRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0incide, line, n, nt, state);
    delete[] line;
    return state;
}
dd::Edge RegisterFactory::ModuloNAdderDebug(vector<int> num1, vector<int> num2){
    assert(num1.size() == num2.size());
    int n = static_cast<int> (num1.size());
    short* line = new short[n];
    gg.lineClear(line, n);
    int nt = n * 3 + 1;//first number, second number, carry and overflow carry.
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    std::function<void (int, int)> lambdaInitState = ExtractedStateInitializer(line, nt, state);
    auto c0indice = [](int i){ return 3 * i;};
    auto a0indice = [](int i){ return 3 * i + 1;};
    auto b0indice = [](int i){ return 3 * i + 2;};
    auto c1indice = [](int i){ return 3 * (i + 1);};
    for(int i = 0; i < n; ++i){
        lambdaInitState(num1[i], a0indice(i));//put first number in register
        lambdaInitState(num2[i], b0indice(i));//put second number in register
    }
    state = RippleAdderGeneral(state, n);
    delete[] line;
    return state;
}
dd::Edge RegisterFactory::ModuloNAdderHalfClassicDebug(lli cnum, vector<int> qnum){
    n = static_cast<int>(base2N.size());
    assert(n == qnum.size());//assumed quantum number represented in n base2 digits.
    short* line = new short[n];
    gg.lineClear(line, n);
    vector<bool> a0clb2 = shor::base2rep(cnum, n);//classical value in base 2
    int nt = 2 * n + 2;//quantum number (n), carries(n), overflow carries(2).
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    std::function<void (int, int)> lambdaInitState = ExtractedStateInitializer(line, nt, state);
    //lambdas used for converting digit (of input numbber) index to qubit index
    //a0 line does not exist. replaced by classic register.
    auto c0indice = [nt](int i){ return nt - 2 - 2 * i;};
    auto a0cnumindice = [a0clb2](int i){return a0clb2.empty() ? 0 : a0clb2[i]; };//return cnum in base 2.
    auto a0Nindice = [Nrep = this->base2N](int i){return Nrep.empty() ? 0 : Nrep[i]; };//return cnum in base 2.
    auto b0indice = [nt](int i){ return nt - 2 - (2 * i + 1);};
    auto c1indice = [nt](int i){ return nt - 2 - 2 * (i + 1);};
    auto tindice = [nt](){return nt - 1;};//second overflow qubit(t in fig 2 of paper)
    for(int i = 0; i < n; ++i){
        lambdaInitState(qnum[i], b0indice(i));//put quantum number in register
    }
    ExtractedRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0cnumindice, line, n, nt, state);
    ExtractedRippleSubtractorHalfClassic(b0indice, c0indice, c1indice, a0Nindice, line, nt, state);
    gg.NotGenOrApply(line, c1indice(n - 1)/*nt - 2*/, nt, &state);
    gg.CNotGenOrApply(line, tindice(), c1indice(n - 1)/*nt - 2*/, nt, &state);
    gg.NotGenOrApply(line, c1indice(n - 1)/*nt - 2*/, nt, &state);
    ExtractedControlledRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0Nindice, tindice(), line, n, nt, state);
    ExtractedRippleSubtractorHalfClassic(b0indice, c0indice, c1indice, a0cnumindice, line, nt, state);
    gg.CNotGenOrApply(line, tindice(), c1indice(n - 1), nt, &state);
    ExtractedRippleAdderHalfClassic(b0indice, c0indice, c1indice, a0cnumindice, line, n, nt, state);
    
    delete[] line;
    return state;
}
void RegisterFactory::ExtractedRippleSubtractorHalfClassic(const std::function<int (int)> &b0indice, const std::function<int (int)> &c0indice, const std::function<int (int)> &c1indice, const std::function<int (int)> &a0indice, short *line, int nt, dd::Edge &state) {
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

/// Used as subtractor in adder modulo N. In this version I am assuming 'a' and 'b' are both qnumbers. Fig 4 of paper.
/// @param cnum classical value to be subtracted
/// @param state state representing qubits
/// @param awhc adder was half classic (in fig 4 of paper)
dd::Edge RegisterFactory::RippleSubtractorHalfClassicGeneral(lli cnum, dd::Edge state, bool awhc){
    vector<bool> a0clb2 = shor::base2rep(cnum, n);//classical value in base 2
    assert(a0clb2.size() == n);
    int nt = awhc ? n * 2 + 1 : n * 3 + 1;//quantum register, carries and overflow qubit
    auto c0indice = [nt, awhc](int i){ return awhc ? nt - 1 - (2 * i) : nt - 1 - (3 * i);};
    auto a0indice = [nt, awhc, a0clb2](int i){ return awhc ? (a0clb2.empty()? 0 : a0clb2[i]) : nt - 1 - (3 * i);};
    auto b0indice = [nt, awhc](int i){ return awhc ? nt - 1 - (2 * i + 1) : nt - 1 - (3 * i + 2);};
    auto c1indice = [nt, awhc](int i){ return awhc ? nt - 1 - 2 * (i + 1) : nt - 1 - 3 * (i + 1);};
    short* line = new short[nt];
    gg.lineClear(line, nt);
    ExtractedRippleSubtractorHalfClassic(a0indice, b0indice, c0indice, c1indice, line, nt, state);
    delete[] line;
    return state;
}
RegisterFactory::~RegisterFactory(){
    delete[] line;
}
