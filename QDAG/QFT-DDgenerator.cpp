//
//  QFT-DDgenerator.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-02.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "IIC-JKU/DDpackage.h"
#include "IIC-JKU/DDcomplex.h"
#include "util.h"
#include "QFT-DDgenerator.hpp"
#include "commonheaders.h"
/*---BEGIN STATE GENERATOR DEFINITIONS---*/
StateGenerator::StateGenerator(dd::Package* dd) {
    this->dd = dd;
}

/// Uniform state generator (1/sqrt(2)^n) |1,1,1,...,1>
/// @param n number of qubits
dd::Edge StateGenerator::dd_UniformState(int n) {
    dd::Edge f = dd::Package::DDone;
    dd::Edge edges[4];
    edges[1] = edges[3] = dd::Package::DDzero;
    
    for (short p = 0; p < n; ++p) {
        edges[0] = f;
        edges[2] = f;
        f = dd->makeNonterminal(p, edges);
    }
    dd::ComplexValue c{1/std::sqrt(pow(2,n)), 0.0 };
    dd::Complex cx = dd->cn.getCachedComplex(c.r,c.i);
    f.w =dd->cn.mulCached(f.w, cx);
    return f;
}
/// Generate 'random' state from square root of 3. digit even:0, digit odd: 1.
/// @param n number of qubits
dd::Edge StateGenerator::dd_Sqrt3State(int n){
    int countodd = 0;
    dd::Edge e_state {nullptr, {nullptr,nullptr}};
    string sqrt3frac = "732050807568877193176604123436845839024";
    for (int i = 0; i < n; ++i){
        if((stoi(std::to_string(sqrt3frac[i])) % 2)){
            ++countodd;
            if(e_state.p == nullptr){
                e_state = dd->makeBasisState(n, i);
                continue;
            }
            e_state = dd->add(e_state, dd->makeBasisState(n, i));
        }
    }
    dd::ComplexValue c{1/sqrt(countodd), 0.0 };//normalize the state.
    dd::Complex cx = dd->cn.getTempCachedComplex(c.r, c.i);
    dd->cn.mul(e_state.w, e_state.w, cx);
    return e_state;
}

/// Generate Random State.
/// @param n num qubits
dd::Edge StateGenerator::dd_RandomState(int n, int seed){
    std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<long> dis0(0.0, pow(2,n) - 1);
    std::uniform_int_distribution<int> dis(0.0, 1.0);
    long r0 = dis0(gen);//at least one of the basis has to have non-zero coefficient.
    int nonzerocount = 0;
    dd::Edge f = dd::Package::DDone;
    dd::Edge edges[4];
    edges[1] = edges[3] = dd::Package::DDzero;
    for (short p = 0; p < n; ++p) {
        bool r = dis(gen);
        bool r0masked = ((r0 >> p) & 1);
        if(r0masked & r){
            edges[2] = f;
            edges[0] = dd::Package::DDzero;
        }
        else if(!r0masked & !r){
            edges[0] = f;
            edges[2] = dd::Package::DDzero;
        }
        else{
            ++nonzerocount;
            edges[2] = f;
            edges[0] = f;
        }
        f = dd->makeNonterminal(p, edges);
    }
    dd::ComplexValue c{1/std::sqrt(pow(2,nonzerocount)), 0.0 };//normalize the state.
    dd::Complex cx = dd->cn.getTempCachedComplex(c.r,c.i);
    dd->cn.mul(f.w, f.w, cx);
    return f;
}

/// return base state. put here to have a central state generator. it just returns the DDpackage base vector.
/// @param n num qubits
/// @param i base index. 0 to 2^n - 1;
dd::Edge StateGenerator::dd_BaseState(int n, lli i) {
    // assert( i >= 0 & i < pow(2,n));
    return dd->makeBasisState(n, i);
}

/// Custom state builder. Function normalizes the final state.
/// @param v complex amplitudes of qubits
/// @param n number of qubits(added as safety feature).
dd::Edge StateGenerator::dd_CustomState(vector<dd::ComplexValue> v, int n) {
    assert (n == log2(static_cast<int>(v.size())));
     dd::Edge state {nullptr, {nullptr, nullptr}};//skeleton to start with.
    dd::Edge temp;
    for (int i = 0; i < v.size(); ++i){
        if(dd::operator != (v[i], {0, 0})){
            dd::Complex vx = dd->cn.getTempCachedComplex(v[i].r,v[i].i);
            temp = dd->makeBasisState(n, i);
            temp.w = dd->cn.mulCached(temp.w, vx);
            //  dd->cn.mul(temp.w, temp.w, vx);TODO: this corrupts data (if used instead of temp.w = dd->cn.mulCached(temp.w, vx)). Understand why.
            if(!state.p){
                state = temp;
            }
            else{
                state = dd->add(temp, state);
            }
        }
    }
    //normalize the state:
    fp summag2 = 0;
    std::for_each(v.begin(), v.end(), [&summag2](dd::ComplexValue x){ summag2 = summag2 +  pow(x.i,2) + pow(x.r,2);});
    dd::Complex tempc = dd->cn.getTempCachedComplex(1/sqrt(summag2), 0);
    dd->cn.mul(state.w, state.w, tempc);
    return state;
}


/*---END STATE GENERATOR DEFINITIONS---*/

/*---BEGIN GATE GENERATOR DEFINITIONS---*/

/// set 'line' for controlled gates.//NOTE: Just for base 2
/// @param line linearray
/// @param t target index
/// @param c0 control index
/// @param c1 control index
/// @param c0p true means c0 is positive control
/// @param c1p true means c1 is positive control
void GateGenerator::lineSet(short *line, int t, int c0, int c1, bool c0p, bool c1p) {
    if(t >= 0)
        line[t] = 2;
    if(c0 >= 0){
        if(c0p)
            line[c0] = 1;
        else
            line[c0] = 0;
    }
    if(c1 >= 0){
        if(c1p)
            line[c1] = 1;
        else
            line[c1] = 0;
    }
}
void GateGenerator::lineSet(short* line, int t, const map<int,bool> &m){
    if(t >= 0)
        line[t] = 2;
    for(auto it = m.cbegin(); it != m.cend(); ++it){
        it->second ? line[it->first] = 1: line[it->first] = 0;
    }
}
/// reset 'line' for controlled gates.//NOTE: Just for base 2
/// @param line linearray
/// @param t target index
/// @param c0 control index
/// @param c1 control index
void GateGenerator::lineReset(short *line, int t, int c0, int c1) {
    if(t >= 0)
        line[t] = -1;
    if(c0 >= 0)
        line[c0] = -1;
    if(c1 >= 0)
        line[c1] = -1;
}
/// reset 'line' for controlled gates, arbitrary number of controls.//NOTE: Just for base 2
/// @param line linearray
/// @param t target index
/// @param m map containing control line indices (first) and whether it is positive or negative (second)
void GateGenerator::lineReset(short* line, int t, const map<int,bool> &m){
    line[t] = -1;
    for(auto it = m.cbegin(); it != m.cend(); ++it){
        line[it->first] = -1;
    }
};
void GateGenerator::lineClear(short *line, int n){
    for(int i = 0; i < n; ++i){
        line[i] = -1;
    }
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
    dd::ComplexValue c{0.5, 0.0 };//used to normalize the state.
    dd::Complex cx = dd->cn.getTempCachedComplex(c.r,c.i);
    ret.w = dd->cn.mulCached(ret.w, cx);
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
/// /inverse rotation gate
/// @param m rotation operator on target qubit
/// @param k  rotation order
void GateGenerator::RInvmatGenerator(dd::Matrix2x2 &m, int k) {
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
        dd::Edge smat = Smatv1(n, i, n - i - 1);
        e_swap = dd->multiply(smat, e_swap);
    }
    return e_swap;
}

/// apply swap gates, one by one to the input state, used in     RegisterFactory::ExponentiateOutReg(). Refer to figure 6 of PRA S1050-2947(96)05707-1.
/// @param nt total number of variables
/// @param v1 vector containing first group of qubits for swap
/// @param v2 vector containing second group of qubits for swap
/// @param state state for operator to be appiled to.
dd::Edge GateGenerator::swapRegistersOnState(int nt, vector<int> v1, vector<int> v2, dd::Edge state) {
    assert (v1.size() == v2.size());
    for(int i = 0; i < v1.size(); ++i){
        state = dd->multiply(Smatv1(nt, v1[i], v2[i]), state);
    }
    return state;
}

/// permute operator, apply swap gates, one by one to the input state
/// @param n number of qubits
/// @param state input state root.
dd::Edge GateGenerator::permuteOperatorOnState(int n, dd::Edge state){//TODO: pass state by value may cause wasted memory. investigate.
    for(int i = 0; i < n/2; ++i){//apply n/2 swap gates.
        state = dd->multiply(Smatv1(n, i, n - i - 1), state);
    }
    return state;
};

/// Controlled-k Not gate. Not gate has arbitrary number of control qubits.
/// @param line gate generator line
/// @param t target index
/// @param nt max index(exclusive)
/// @param m map containing control line indices (first) and whether it is positive or negative (second). map chosen over hash map since num qubits is small.
/// @param state if provided, gate will be applied to this state
/// @return Controlled-k Not gate.
dd::Edge GateGenerator::CKNotGenOrApply(short* line, int t, const map<int,bool> &m, int nt,  dd::Edge* state){
    lineSet(line, t, m);
    dd::Edge e = dd->makeGateDD(Xmat, nt, line);
    lineReset(line, t, m);
    if(state){
        *state = dd->multiply(e, *state);
    }
    return e;
}
/// Create Toffoli gate if not computed already.
/// @param line gate generator line
/// @param c0 first control index
/// @param c1 second control index
/// @param t target index
/// @param nt max index(exclusive)
/// @param state optional state value for gate to be applied to before return.
/// @param c0p true means c0 is positive control
/// @param c1p true means c1 is positive control
dd::Edge GateGenerator::ToffoliGenOrApply(short* line, int t, int c0, int c1, int nt, dd::Edge* state, bool c0p, bool c1p){
    lineSet(line, t, c0, c1, c0p, c1p);
    dd::Edge e = dd->makeGateDD(Xmat, nt, line);
    lineReset(line, t, c0, c1);
    if(state){
        *state = dd->multiply(e, *state);
    }
    return e;
}

/// Create C-Not gate//Can be combined with ToffoliGenerator(), left for code readability.
/// @param line gate generator line
/// @param t target index
/// @param c control index
/// @param nt max indix(exclusive)
/// @param state optional state value for gate to be applied to before return.
/// @param cp true means control is positive
dd::Edge GateGenerator::CNotGenOrApply(short* line, int t, int c, int nt, dd::Edge* state,  bool cp){
    lineSet(line, t, c, -1, cp);
    dd::Edge e = dd->makeGateDD(Xmat, nt, line);
    lineReset(line, t, c);
    if(state){
        *state = dd->multiply(e, *state);
    }
    return e;
}

/// Create Not-gate//Can be combined with CNotGenerator() or ToffoliGenerator(), left for code readability.
/// @param line gate generator line
/// @param t target index
/// @param nt max index(exclusive)
/// @param state optional state value for gate to be applied to before return.
dd::Edge GateGenerator::NotGenOrApply(short* line, int t, int nt, dd::Edge* state){
    lineSet(line, t);
    dd::Edge e = dd->makeGateDD(Xmat, nt, line);
    lineReset(line, t);
    if(state){
        *state = dd->multiply(e, *state);
    }
    return e;
}
/// Create Hadamard. be combined with CNotGenerator() or ToffoliGenerator(), left for code readability.
/// @param line gate generator line
/// @param t target index
/// @param nt max index(exclusive)
/// @param state optional state value for gate to be applied to before return.
dd::Edge GateGenerator::HadGenOrApply(short *line, int t, int nt, dd::Edge* state){
    lineSet(line, t);
    dd::Edge e = dd->makeGateDD(Hmat, nt, line);
    lineReset(line, t);
    if(state){
        *state = dd->multiply(e, *state);
    }
    return e;
}
/*---END GATE GENERATOR DEFINITIONS---*/
