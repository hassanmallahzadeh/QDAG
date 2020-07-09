//
//  QFT.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-06-09.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "util.h"
#include "IIC-JKU/DDcomplex.h"
#include "QFT-DDgenerator.hpp"
#include "QFT.hpp"
using namespace std::chrono;
using std::cout;
using std::endl;
QFT::QFT(dd::Package* dd){
    this->dd = dd;
}
/// QFT, make the overal operator (interminnet 'line' gates constructed, applied to answer incrementally)
///  then apply on input state.
/// uncomment export2Dot and cout lines to monitor bdd structure evolution. comment to disable
/// @param n number of qubits
/// @param perm where and if apply permutation gate.
dd::Edge QFT::dd_QFTV1(int n, PERM_POS perm) {
    dd::Edge e_ans = dd->makeIdent(0, n-1);//HM: this can be taken out, here for clarity
    GateGenerator  gg = GateGenerator(dd);
    
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    
    for (int i = 0; i < n; ++i){//HM: for each of n qubits do:
        gg.lineSet(line,i,-1);
        dd::Edge e_line = dd->makeGateDD(Hmat, n, line);//HM: start by hadamard as in circuit
        dd::Matrix2x2 Rmat;//HM: put rotation gates in place
        int j = 0;
        for (; j < n - i - 1; ++j){
            {
                if(m_ord == REG_C_T)
                    gg.lineSet(line, i, i + j + 1);
                else
                    gg.lineSet(line, i + j + 1, i);
            }
            gg.RmatGenerator(Rmat, j+2);
            dd::Edge e_temp = dd->makeGateDD(Rmat, n, line);
            e_line = dd->multiply(e_temp, e_line);
            dd->export2Dot(e_line, "qft-l-"+std::to_string(i)+"-"+std::to_string(j)+".dot",false, true);
            {
                if(m_ord == REG_C_T)
                    gg.lineReset(line, i, i + j + 1);
                else
                    gg.lineReset(line, i + j + 1, i);
            }
            
        }
        e_ans = dd->multiply(e_line, e_ans);//HM: apply current line to final operator
        gg.lineReset(line,i,n - i - 1);//HM: reset line array
    }
    if(perm == BEG_PERM){
        e_ans = dd->multiply(e_ans, gg.permuteOperator(n));
    }
    else if(perm == END_PERN){
        e_ans = dd->multiply(gg.permuteOperator(n),e_ans);//this is equivanent to reversing variable order.
    }
    delete[] line;
    return e_ans;
}

/// QFT, apply the gates one by one to the input state. avoid making the overal operator first.
/// @param n number of qubits
/// @param state input state
/// @param perm where and if apply permutation gate.
dd::Edge QFT::dd_QFTV2(int n, dd::Edge state, PERM_POS perm) {
    GateGenerator  gg = GateGenerator(dd);
    dd::Edge e_ans = state;
    if(perm == BEG_PERM){
        e_ans = dd->multiply(gg.permuteOperator(n),e_ans);
    }
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    for (int i = 0; i < n; ++i){//HM: for each of n qubits do:
        gg.lineSet(line,i,-1);
        dd::Edge e_line = dd->makeGateDD(Hmat, n, line);//HM: start by hadamard as in circuit
        e_ans = dd->multiply(e_line ,e_ans);
        dd::Matrix2x2 Rmat;//HM: put rotation gates in place
        for (int j = 0; j < n - i - 1; ++j){
            if(m_ord == REG_C_T)
                gg.lineSet(line, i, i + j + 1);
            else
                gg.lineSet(line, i + j + 1, i);
            gg.RmatGenerator(Rmat, j+2);
            e_ans = dd->multiply(dd->makeGateDD(Rmat, n, line), e_ans);
            if(m_ord == REG_C_T)
                gg.lineReset(line, i, i + j + 1);
            else
                gg.lineReset(line, i + j + 1, i);
        }
        gg.lineReset(line,i,n - i - 1);//HM: reset line array
    }
    if(perm == END_PERN)
        e_ans = dd->multiply(gg.permuteOperator(n),e_ans);
    
    delete[] line;
    return e_ans;
}
/// QFT, make the overal operator, then apply on input state. (interminnet 'line' gates not constructed, multiply all gates in order of appearing)
/// uncomment export2Dot and cout lines to monitor bdd structure evolution. comment to disable
/// @param n number of qubits
/// @param perm where and if apply permutation gate.
dd::Edge QFT::dd_QFTV3(int n , PERM_POS perm){
    dd::Edge e_ans = dd->makeIdent(0, n-1);//HM: this can be taken out, here for clarity
    GateGenerator  gg = GateGenerator(dd);
    if(perm == BEG_PERM){
        e_ans = dd->multiply(e_ans, gg.permuteOperator(n));
    }
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    
    for (int i = 0; i < n; ++i){//HM: for each of n qubits do:
        gg.lineSet(line,i,-1);
        e_ans = dd->multiply(dd->makeGateDD(Hmat, n, line), e_ans);//HM: start by hadamard as in circuit
        dd::Matrix2x2 Rmat;//HM: put rotation gates in place
        //        dd->export2Dot(e_ans, "qft-h-"+std::to_string(i)+".dot",false, true);
        
        int j = 0;
        for (; j < n - i - 1; ++j){
            if(m_ord == REG_C_T)
                gg.lineSet(line, i, i + j + 1);
            else
                gg.lineSet(line, i + j + 1, i);
            gg.RmatGenerator(Rmat, j+2);
            dd::Edge e_rot = dd->makeGateDD(Rmat, n, line);
            e_ans = dd->multiply(e_rot, e_ans);
            if(m_ord == REG_C_T)
                gg.lineReset(line, i, i + j + 1);
            else
                gg.lineReset(line, i + j + 1, i);
            
        }
        gg.lineReset(line,i,n - i - 1);//HM: reset line array
    }
    if(perm == END_PERN){
        e_ans = dd->multiply(gg.permuteOperator(n),e_ans);//this is equivanent to reversing variable order.
    }
    delete[] line;
    return e_ans;
}
