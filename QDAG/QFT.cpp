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
#include "QFT-Measurement.hpp"
using namespace std::chrono;
using std::cout;
using std::endl;

const int QFT::posControl = 1;//apply rotation gates on stage i if qubit measrement on stage i had result posControl. refer to Mermin, fig 3.3. Would be 0, 1 or 2 for base 3, etc.

QFT::QFT(dd::Package* dd){
    assert(dd::RADIX == 2);//parts of this code are not properly made for base != 2.
    this->dd = dd;
}
/// QFT, make the overal operator (interminnet 'line' gates constructed, applied to answer incrementally)
///  then apply on input state.
/// uncomment export2Dot and cout lines to monitor bdd structure evolution. comment to disable
/// @param n number of qubits
/// @param perm where and if apply permutation gate.
dd::Edge QFT::dd_QFTV1(int n, PERM_POS perm) {
    dd::Edge e_ans = dd->makeIdent(0, n-1);// this can be taken out, here for clarity
    GateGenerator  gg = GateGenerator(dd);
    
    short *line = new short[n]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    
    for (int i = 0; i < n; ++i){// for each of n qubits do:
        gg.lineSet(line,i,-1);
        dd::Edge e_line = dd->makeGateDD(Hmat, n, line);// start by hadamard as in circuit
        dd::Matrix2x2 Rmat;// put rotation gates in place
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
        e_ans = dd->multiply(e_line, e_ans);// apply current line to final operator
        gg.lineReset(line,i,n - i - 1);// reset line array
    }
    if(perm == BEG_PERM){
        e_ans = gg.permuteOperatorOnState(n, e_ans);
    }
    else if(perm == END_PERM){
        e_ans = gg.permuteOperatorOnState(n, e_ans);//this is equivanent to reversing variable order.
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
        e_ans = gg.permuteOperatorOnState(n, e_ans);
    }
    else if(perm == END_PERM)
        cout<<"Warning: permutation gate has to be applied at beginning for correct result.\n";
    short *line = new short[n]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    for (int i = 0; i < n; ++i){// for each of n qubits do:
        gg.lineSet(line,i,-1);
        dd::Edge e_line = dd->makeGateDD(Hmat, n, line);// start by hadamard as in circuit
        e_ans = dd->multiply(e_line ,e_ans);
        dd::Matrix2x2 Rmat;// put rotation gates in place
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
    }
    if(perm == END_PERM){
        e_ans = gg.permuteOperatorOnState(n, e_ans);
    }
    
    delete[] line;
    return e_ans;
}
/// QFT, make the overal operator, then apply on input state. (interminnet 'line' gates not constructed, multiply all gates in order of appearing)
/// uncomment export2Dot and cout lines to monitor bdd structure evolution. comment to disable
/// @return overall QFT decition diagram.
/// @param n number of qubits
/// @param perm where and if apply permutation gate.
dd::Edge QFT::dd_QFTV3(int n , PERM_POS perm){
    dd::Edge e_ans = dd->makeIdent(0, n-1);// this can be taken out, here for clarity
    GateGenerator  gg = GateGenerator(dd);
    if(perm == BEG_PERM){
        e_ans = gg.permuteOperatorOnState(n, e_ans);
    }
    short *line = new short[n]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    
    for (int i = 0; i < n; ++i){// for each of n qubits do:
        gg.lineSet(line,i,-1);
        e_ans = dd->multiply(dd->makeGateDD(Hmat, n, line), e_ans);// start by hadamard as in circuit
        dd::Matrix2x2 Rmat;// put rotation gates in place
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
        gg.lineReset(line,i,n - i - 1);// reset line array
    }
    if(perm == END_PERM){
        e_ans = gg.permuteOperatorOnState(n, e_ans);//this is equivanent to reversing variable order.
    }
    delete[] line;
    return e_ans;
}

/// QFT by permuting rows and columns of operator matrix to be able to apply permutation gate at the end.
/// @return state after act of QFT.
/// @param n number of qubits
/// @param state input state root edge
/// @param perm where and if apply the permutation operator. is here for consistency with other versions of the function
dd::Edge QFT::dd_QFTV4(int n, dd::Edge state, PERM_POS perm) {
    
    GateGenerator  gg = GateGenerator(dd);
    dd::Edge e_ans = state;
    
    if(perm == BEG_PERM){
        e_ans = gg.permuteOperatorOnState(n, e_ans);
        cout<<"Warning: permutation gate has to be applied at end for correct result.\n";
    }
    short *line = new short[n]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    for (int i = n - 1; i >= 0; --i){// for each of n qubits do:
        gg.lineSet(line, i, -1);
        e_ans = dd->multiply(dd->makeGateDD(Hmat, n, line) ,e_ans); //start by hadamard as in circuit
        
        dd::Matrix2x2 Rmat;// put rotation gates in place
        for (int j = i - 1; j >= 0; --j){
            if(m_ord == REG_C_T)
                gg.lineSet(line, i, j);
            else
                gg.lineSet(line, j, i);
            
            gg.RmatGenerator(Rmat, i - j + 1);
            e_ans = dd->multiply(dd->makeGateDD(Rmat, n, line), e_ans);
            if(m_ord == REG_C_T)
                gg.lineReset(line, i, j);
            else
                gg.lineReset(line, j, i);
        }
    }
    if(perm == END_PERM)
        e_ans = gg.permuteOperatorOnState(n, e_ans);
    
    delete[] line;
    return e_ans;
}



/// QFT with Griffiths-Niu scheme.
/// @return state after act of QFT.
/// @param n number of qubits
/// @param state input state root edge
/// @param perm where and if apply the permutation operator.
dd::Edge QFT::dd_QFTGNV1(int n, dd::Edge state, PERM_POS perm, engine& unrg){
//assert(m_ord == INV_C_T);

    GateGenerator  gg = GateGenerator(dd);
    dd::Edge e_ans = state;
    Measurement mm = Measurement(dd);
    if(perm == BEG_PERM){
        e_ans = gg.permuteOperatorOnState(n, e_ans);
dd->garbageCollect();//without garbage collect crash happens on more than 17 qubits on my computer.//TODO: go after finding why crash happens.
    }
    
    short *line = new short[n]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    for (int i = 0; i < n; ++i){// for each of n qubits do:
        gg.lineSet(line, i, -1);
        e_ans = dd->multiply(dd->makeGateDD(Hmat, n, line) ,e_ans); //start by hadamard as in circuit
        gg.lineReset(line, i, -1);
        int res = mm.Measure(e_ans, n, i, unrg);
        if(res == posControl){
            dd::Matrix2x2 Rmat;// put rotation gates in place
            for (int j = 0; j < n - i - 1; ++j){
                //activate rotation gates
                gg.lineSet(line, i + j + 1, -1);
                gg.RmatGenerator(Rmat, j+2);
                e_ans = dd->multiply(dd->makeGateDD(Rmat, n, line), e_ans);
                gg.lineReset(line, i + j + 1, -1);
            }
            gg.lineReset(line, i, n - i - 1);// likely redundant. reset line array
        }
    }
    if(perm == END_PERM){
        dd->garbageCollect();
      e_ans = gg.permuteOperatorOnState(n, e_ans);
    cout<<"Warning: permutation gate has to be applied at the beginning for correct result.\n";
    }

    delete[] line;
    return e_ans;

}
/// QFT with Griffiths-Niu scheme. permuted QFT operator.
/// @return state after act of QFT.
/// @param n number of qubits
/// @param state input state root edge
dd::Edge QFT::dd_QFTGNV2(int n, dd::Edge state, PERM_POS perm ,engine& unrg) {
    assert(m_ord == INV_C_T);//to be taken out

    GateGenerator  gg = GateGenerator(dd);
    dd::Edge e_ans = state;
    if(perm == BEG_PERM){
          e_ans = gg.permuteOperatorOnState(n , e_ans);
          cout<<"Warning: permutation gate has to be applied at end for correct result.\n";
      }
    Measurement mm = Measurement(dd);
    short *line = new short[n]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    for (int i = n - 1; i >= 0; --i){// for each of n qubits do:
        gg.lineSet(line, i, -1);
        e_ans = dd->multiply(dd->makeGateDD(Hmat, n, line) ,e_ans); //start by hadamard as in circuit
        gg.lineReset(line, i, -1);
        int res = mm.Measure(e_ans, n, i, unrg);
        if(res == posControl){
            dd::Matrix2x2 Rmat;// put rotation gates in place
            for (int j = i - 1; j >= 0; --j){
                //activate rotation gates
                gg.lineSet(line, j, -1);
                gg.RmatGenerator(Rmat, i - j + 1);
                e_ans = dd->multiply(dd->makeGateDD(Rmat, n, line), e_ans);
                gg.lineReset(line, j, -1);
            }
          //  gg.lineReset(line, i, -1);// reset line array
        }
    }
    if(perm == END_PERM){
           e_ans = gg.permuteOperatorOnState(n, e_ans);
       }
    delete[] line;
    return e_ans;
}
/// Applies QFT with GN Scheme to qubit indices specified in vector i.
/// @return outcomes measurement outcomes during QFT.
/// @param n number of qubits
/// @param state input state root edge
/// @param perm where and if apply the permutation operator.
/// @param indice indice of qubits for applying QFT
vector<int> QFT::dd_QFTGNV1(int n, dd::Edge &state, PERM_POS perm, engine& unrg, vector<int> indice){
    GateGenerator  gg = GateGenerator(dd);
    if(perm == BEG_PERM){
        state = gg.permuteOperatorOnState(n, state);
//dd->garbageCollect();//without garbage collect crash happens on more than 17 qubits on my computer.//TODO: go after finding why crash happens.
    }
    short *line = new short[n]{};// set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    vector<int> measurments;
    Measurement mm = Measurement(dd);
    for(int i = 0; i < indice.size(); ++i){
          gg.lineSet(line, indice[i]);
        state = dd->multiply(dd->makeGateDD(Hmat, n, line) ,state); //start by hadamard as in circuit
        gg.lineReset(line, indice[i]);
        int res = mm.Measure(state, n, indice[i], unrg);
        measurments.push_back(res);
        if(res == posControl){
            dd::Matrix2x2 Rmat;// put rotation gates in place
            for (int j = 0; j < indice.size() - i - 1; ++j){
                gg.lineSet(line, indice[i + j + 1], -1);
                gg.RmatGenerator(Rmat, j+2);
                state = dd->multiply(dd->makeGateDD(Rmat, n, line), state);
                gg.lineReset(line, indice[i + j + 1], -1);
            }
        }
    }
    if(perm == END_PERM){
    //    dd->garbageCollect();
      state = gg.permuteOperatorOnState(n, state);
    cout<<"Warning: permutation gate has to be applied at the beginning for correct result.\n";
    }
 //   dd->garbageCollect();
    return measurments;
}
