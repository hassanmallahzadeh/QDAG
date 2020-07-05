//
//  QFT-DDgenerator.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-02.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "DDpackage.h"
#include "DDcomplex.h"
#include "util.h"
#include "QFT-DDgenerator.hpp"
using CN = dd::ComplexNumbers;
StateGenerator::StateGenerator(dd::Package* dd) {
    this->dd = dd;
}


/// Uniform state generator (1/sqrt(2)^n) |1,1,1,...,1>
/// @param n number of qubits
dd::Edge StateGenerator::dd_UniformState(int n) {
    dd::Edge e_state = dd->makeBasisState(n, 0);
    for (int i = 1; i < pow(2,n); ++i){
        e_state = dd->add(e_state, dd->makeBasisState(n, i));
    }
    dd::ComplexValue c{1/std::sqrt(pow(2,n)), 0.0 };
    dd::Complex cx = dd->cn.getTempCachedComplex(c.r,c.i);
    e_state.w =dd->cn.mulCached(e_state.w, cx);
    return e_state;
}
/// Generate 'random' state from square root of 3. digit even:0, digit odd: 1.
/// @param n number of qubits
dd::Edge StateGenerator::dd_Sqrt3State(int n){//TODO: fix: function breaks after 9 steps (return negative for integer).
    fp sqrt3 = sqrt(3);
    dd::Edge e_state = dd->makeBasisState(n, 0);//3 gives 1.
    sqrt3 = sqrt3 * 10;
    for (int i = 1; i < pow(2,n); ++i){
        if((int)sqrt3 % 2)//bit is 1, so add
            e_state = dd->add(e_state, dd->makeBasisState(n, i));
        sqrt3 = sqrt3 * 10;
    }
    dd::ComplexValue c{1/std::sqrt(pow(2,n)), 0.0 };//normalize the state.
    dd::Complex cx = dd->cn.getTempCachedComplex(c.r,c.i);
    e_state.w =dd->cn.mulCached(e_state.w, cx);
    return e_state;
}

/// Generate Random State. TODO: replace with c++11 random number generator.
/// @param n num qubits
dd::Edge StateGenerator::dd_RandomState(int n){
    dd::Edge e_state = dd->makeBasisState(n, 0);//assume first one is 1 (TODO: fix)
    srand (static_cast<unsigned int>(time(nullptr)));
    for (int i = 1; i < pow(2,n); ++i){
        if(rand() % 2)//bit is 1, so add
            e_state = dd->add(e_state, dd->makeBasisState(n, i));
    }
    dd::ComplexValue c{1/std::sqrt(pow(2,n)), 0.0 };//normalize the state.
    dd::Complex cx = dd->cn.getTempCachedComplex(c.r,c.i);
    e_state.w =dd->cn.mulCached(e_state.w, cx);
    return e_state;
}

/// return base state. put here to have a central state generator. it just returns the DDpackage base vector.
/// @param n num qubits
/// @param i base index. 0 to 2^n - 1;
dd::Edge StateGenerator::dd_BaseState(int n, int i) {
    assert( i >= 0 & i < pow(2,n));
    return dd->makeBasisState(n, i);
}
/// set 'line' for controlled gates.
/// @param line linearray
/// @param t target index
/// @param c control index
void GateGenerator::lineSet(short *line, int t, int c) { 
    if(t >= 0)
        line[t] = 2;
    if(c >= 0)
        line[c] = 1;
}
/// reset 'line' for controlled gates.
/// @param line linearray
/// @param t target index
/// @param c control index
void GateGenerator::lineReset(short *line, int t, int c) { 
    if(t >= 0)
        line[t] = -1;
    if(c >= 0)
        line[c] = -1;;
}

GateGenerator::GateGenerator(dd::Package *dd) { 
    this->dd = dd;
}
/// Smatv1: Swap matrix construct using Pauli X, Y, Z, I matrices.
/// @param n number of bits
/// @param b1 first bit in swap
/// @param b2 second bit in swap
dd::Edge GateGenerator::Smatv1(int n, int b1, int b2) { 
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    line[b1] = 2;
    dd::Edge x1_gate = dd->makeGateDD(Xmat, n, line);
    dd::Edge y1_gate = dd->makeGateDD(Ymat, n, line);
    dd::Edge z1_gate = dd->makeGateDD(Zmat, n, line);
    line[b1] = -1;
    line[b2] = 2;
    dd::Edge x2_gate = dd->makeGateDD(Xmat, n, line);
    dd::Edge y2_gate = dd->makeGateDD(Ymat, n, line);
    dd::Edge z2_gate = dd->makeGateDD(Zmat, n, line);
    dd::Edge iMatnxn = dd->makeIdent(0, n-1);
    dd::Edge x_gate = dd->multiply(x1_gate, x2_gate);//make 2 qubit dd
    dd::Edge y_gate = dd->multiply(y1_gate, y2_gate);//make 2 qubit dd
    dd::Edge z_gate = dd->multiply(z1_gate, z2_gate);//make 2 qubit dd
    dd::Edge temp1 = dd->add(iMatnxn, x_gate);
    dd::Edge temp2 = dd->add(y_gate, z_gate);
    dd::Edge ret = dd->add(temp1, temp2);
    // ret = ret * 2; HM: TODO: code has to be added for this.
    delete[] line;
    return ret;;
}
/// Smatv2: Swap matrix construct using C-Not gates.
/// @param n number of bits
/// @param b1 first bit in swap
/// @param b2 second bit in swap
dd::Edge GateGenerator::Smatv2(int n, int b1, int b2){
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    line[b1] = 2;//HM: set target and control
    line[b2] = 1;
    dd::Edge c12_gate = dd->makeGateDD(Xmat, n, line);
    line[b1] = 1;
    line[b2] = 2;
    dd::Edge c21_gate = dd->makeGateDD(Xmat, n, line);
    dd::Edge temp_gate = dd->multiply(c21_gate, c12_gate);
    dd::Edge s_gate = dd->multiply(temp_gate, c21_gate);
    delete[] line;
    return s_gate;
}
/// /rotation gate
/// @param m rotation operator on target qubit
/// @param k  rotation order
void GateGenerator::RmatGenerator(dd::Matrix2x2 &m, int k) { 
    m[0][0] = { 1, 0 };
    m[0][1] = { 0, 0 };
    m[1][0] = { 0, 0 };
    fp angle = 2 * dd::PI/pow(2,k);
    m[1][1] = { cos(angle), sin(angle) };;
}
/// permute operator. can be placed at beginning or end of circuit
/// @param n number of variables
dd::Edge GateGenerator::permuteOperator(int n) { 
    dd::Edge e_swap = dd->makeIdent(0, n-1);
    for(int i = 0; i < n/2; ++i){//apply n/2 swap gates.
        e_swap = dd->multiply(Smatv2(n, i, n - i - 1), e_swap);
    }
    return e_swap;;
}
