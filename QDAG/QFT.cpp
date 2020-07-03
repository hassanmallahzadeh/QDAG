//
//  QFT.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-06-09.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "util.h"
#include "DDcomplex.h"
#include "QFT.hpp"
using namespace std::chrono;
using std::cout;
using std::endl;
/// Smatv1: Swap matrix construct using Pauli X, Y, Z, I matrices.
/// @param dd source package instance
/// @param n number of bits
/// @param b1 first bit in swap
/// @param b2 second bit in swap
dd::Edge Smatv1(dd::Package *dd, int n, int b1, int b2){
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
    return ret;
}
/// Smatv2: Swap matrix construct using C-Not gates.
/// @param dd source package instance
/// @param n number of bits
/// @param b1 first bit in swap
/// @param b2 second bit in swap
dd::Edge Smatv2(dd::Package *dd, int n, int b1, int b2){
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
void RmatGenerator(dd::Matrix2x2 &m, int k){
    m[0][0] = { 1, 0 };
    m[0][1] = { 0, 0 };
    m[1][0] = { 0, 0 };
    fp angle = 2 * dd::PI/pow(2,k);
    m[1][1] = { cos(angle), sin(angle) };
}
/// set 'line' for controlled gates.
/// @param line linearray
/// @param t target index
/// @param c control index
void lineSet(short* line, int t,int c){
    if(t >= 0)
        line[t] = 2;
    if(c >= 0)
        line[c] = 1;
}
/// reset 'line' for controlled gates.
/// @param line linearray
/// @param t target index
/// @param c control index
void lineReset(short* line, int t,int c){
    if(t >= 0)
        line[t] = -1;
    if(c >= 0)
        line[c] = -1;
}
/// permute operator. can be placed at beginning or end of circuit
/// @param dd package pointer
/// @param n number of variables
dd::Edge permuteOperator(dd::Package *dd, int n){
    dd::Edge e_swap = dd->makeIdent(0, n-1);
    for(int i = 0; i < n/2; ++i){//apply n/2 swap gates.
        e_swap = dd->multiply(Smatv2(dd, n, i, n - i - 1), e_swap);
    }
    return e_swap;
}
/// QFT, make the overal operator, then apply on input state.
/// uncomment export2Dot and cout lines to monitor bdd structure evolution. comment to disable
/// @param dd package pointer
/// @param n number of qubits
/// @param perm where and if apply permutation gate.
dd::Edge dd_QFTV1(dd::Package *dd, int n,  PERM_POS perm){
    dd::Edge e_ans = dd->makeIdent(0, n-1);//HM: this can be taken out, here for clarity
    
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    std::ofstream datafile;
    datafile.open ("timers.txt");
    time_point<system_clock> start, end;
    for (int i = 0; i < n; ++i){//HM: for each of n qubits do:
        start = system_clock::now();
        lineSet(line,i,-1);
        dd::Edge e_line = dd->makeGateDD(Hmat, n, line);//HM: start by hadamard as in circuit
        dd::Matrix2x2 Rmat;//HM: put rotation gates in place
        int j = 0;
        for (; j < n - i - 1; ++j){
            //  dd->export2Dot(e_line, "qft-l-"+std::to_string(i)+"-"+std::to_string(j)+".dot",false, true);
            lineSet(line, -1, i + j + 1);
            RmatGenerator(Rmat, j+2);
            //     cout<<"gate counter: "<<i<<", "<<j<<", "<<dd->size(e_line)<< endl;
            //    dd::Edge e_temp = dd->makeGateDD(Rmat, n, line);
            //     dd->export2Dot(e_temp, "qft-r-"+std::to_string(i)+"-"+std::to_string(j)+".dot",false, true);
            e_line = dd->multiply(dd->makeGateDD(Rmat, n, line), e_line);
            lineReset(line, -1, i + j + 1);
        }
        // cout<<"bit: "<<i<<", gate: "<< j << ", size: "<<dd->size(e_line) << endl;
        //   dd->export2Dot(e_line, "qft-l-"+std::to_string(i)+".dot",false, true);
        e_ans = dd->multiply(e_line, e_ans);//HM: apply current line to final operator
        //    dd->export2Dot(e_ans, "qft-a-"+std::to_string(i)+".dot",false, true);
        lineReset(line,i,n - i - 1);//HM: reset line array
        end = system_clock::now();
        duration<float> elapsed_seconds = end - start;
        datafile << i << " "<< elapsed_seconds.count() <<endl;
    }
    if(perm == BEG_PERM){
        e_ans = dd->multiply(e_ans, permuteOperator(dd,n));
    }
    else if(perm == END_PERN){
        e_ans = dd->multiply(permuteOperator(dd,n),e_ans);//this is equivanent to reversing variable order.
    }
    //   dd->export2Dot(e_ans, "qft-f.dot",false, true);
    datafile.close();
    delete[] line;
    return e_ans;
}
/// QFT, apply the gates one by one to the input state. avoid making the overal operator first.
/// @param dd package pointer
/// @param n number of qubits
/// @param state input state
/// @param perm where and if apply permutation gate.
dd::Edge dd_QFTV2(dd::Package *dd, int n, dd::Edge state, PERM_POS perm){
    dd::Edge e_ans = state;
    if(perm == BEG_PERM){
        e_ans = dd->multiply(permuteOperator(dd,n),e_ans);
    }
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    for (int i = 0; i < n; ++i){//HM: for each of n qubits do:
        lineSet(line,i,-1);
        dd::Edge e_line = dd->makeGateDD(Hmat, n, line);//HM: start by hadamard as in circuit
        e_ans = dd->multiply(e_line ,e_ans);
        dd::Matrix2x2 Rmat;//HM: put rotation gates in place
        for (int j = 0; j < n - i - 1; ++j){
            lineSet(line, -1, i + j + 1);
            RmatGenerator(Rmat, j+2);
            e_ans = dd->multiply(dd->makeGateDD(Rmat, n, line), e_ans);
            lineReset(line, -1, i + j + 1);
        }
        lineReset(line,i,n - i - 1);//HM: reset line array
    }
    if(perm == END_PERN)
        e_ans = dd->multiply(permuteOperator(dd,n),e_ans);
    
    delete[] line;
    return e_ans;
}
/// QFT in main.
int main(){
    
    if(false){//make 'true' to investigate a single QFT (fixed number of bits)
        auto* dd = new dd::Package;
        int n = 5;
        
        //
        //         dd::Edge state3 = dd->makeBasisState(2, 3);
        //         dd::Edge state0 = dd->makeBasisState(n, 0);
        dd::Edge state1 = dd->makeBasisState(n, 1);
        //        dd::Edge state2 = dd->makeBasisState(n, 2);
        //        dd::Edge state3 = dd->makeBasisState(n, 3);
        
        dd::Edge state = state1;
        cout << "Input vector:\n";
        dd->printVector(state);
        dd->export2Dot(state, "state5before.dot",true, true);
        //  dd->printVector(state);
        //  dd::Edge e_swap =permuteOperator(dd,n);//HM: this can be taken out, here for clarity
        state = dd_QFTV2(dd, n, state,NO_PERM);
        dd->export2Dot(state, "state5after.dot",true, true);
        cout << "Output vector:\n";
        dd->printVector(state);
        delete dd;
    }
    
    if(true){//make 'true' to plot the graph (time and node count) as function of varible number.
        int N = 20;//points for graph
        int offset = 2;//starting num of qubits
        float a_et[2][N];//execution time
        int nodecounter[N];
        for (int i = 0; i < N; ++i){
            //Initialize package
            auto* dd = new dd::Package;
            time_point<system_clock> start, end;
            start = system_clock::now();
            dd::Edge e_state = dd->makeBasisState(offset + i, offset + i - 1);
//            cout << "Input vector "<<i<<" bits:\n";
//            dd->printVector(e_state);
            cout<<dd->size(e_state)<<endl;
            e_state = dd_QFTV2(dd, offset + i, e_state, NO_PERM);
//            cout << "Output vector "<<i<<" bits:\n";
//            dd->printVector(e_state);
            nodecounter[i] = dd->size(e_state);
            end = system_clock::now();
            duration<float> elapsed_seconds = end - start;
            a_et[0][i] = offset+i;
            a_et[1][i] = elapsed_seconds.count();
            delete dd;
        }
        std::ofstream datafile;
        datafile.open ("timers.txt");
        if(datafile.is_open()){
            for (int i = 0; i < N; ++i){
                datafile << a_et[0][i] << " "<< a_et[1][i] <<endl;
            }
            datafile.close();
        }
        datafile.open ("counters.txt");
        if(datafile.is_open()){
            for (int i = 0; i < N; ++i){
                datafile << a_et[0][i] << " "<< nodecounter[i] <<endl;
            }
            datafile.close();
        }
    }
    return 0;
}
