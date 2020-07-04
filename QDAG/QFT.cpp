//
//  QFT.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-06-09.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "util.h"
#include "DDcomplex.h"
#include "QFT-DDgenerator.hpp"
#include "QFT.hpp"
using namespace std::chrono;
using std::cout;
using std::endl;
/// QFT, make the overal operator, then apply on input state.
/// uncomment export2Dot and cout lines to monitor bdd structure evolution. comment to disable
/// @param dd package pointer
/// @param n number of qubits
/// @param perm where and if apply permutation gate.
dd::Edge dd_QFTV1(dd::Package *dd, int n,  PERM_POS perm){
    dd::Edge e_ans = dd->makeIdent(0, n-1);//HM: this can be taken out, here for clarity
    GateGenerator  gg = GateGenerator(dd);

    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    std::ofstream datafile;
    datafile.open ("timers.txt");
    time_point<system_clock> start, end;
    for (int i = 0; i < n; ++i){//HM: for each of n qubits do:
        start = system_clock::now();
        gg.lineSet(line,i,-1);
        dd::Edge e_line = dd->makeGateDD(Hmat, n, line);//HM: start by hadamard as in circuit
        dd::Matrix2x2 Rmat;//HM: put rotation gates in place
        int j = 0;
        for (; j < n - i - 1; ++j){
            //  dd->export2Dot(e_line, "qft-l-"+std::to_string(i)+"-"+std::to_string(j)+".dot",false, true);
            gg.lineSet(line, -1, i + j + 1);
            gg.RmatGenerator(Rmat, j+2);
            //     cout<<"gate counter: "<<i<<", "<<j<<", "<<dd->size(e_line)<< endl;
            //    dd::Edge e_temp = dd->makeGateDD(Rmat, n, line);
            //     dd->export2Dot(e_temp, "qft-r-"+std::to_string(i)+"-"+std::to_string(j)+".dot",false, true);
            e_line = dd->multiply(dd->makeGateDD(Rmat, n, line), e_line);
            gg.lineReset(line, -1, i + j + 1);
        }
        // cout<<"bit: "<<i<<", gate: "<< j << ", size: "<<dd->size(e_line) << endl;
        //   dd->export2Dot(e_line, "qft-l-"+std::to_string(i)+".dot",false, true);
        e_ans = dd->multiply(e_line, e_ans);//HM: apply current line to final operator
        //    dd->export2Dot(e_ans, "qft-a-"+std::to_string(i)+".dot",false, true);
        gg.lineReset(line,i,n - i - 1);//HM: reset line array
        end = system_clock::now();
        duration<float> elapsed_seconds = end - start;
        datafile << i << " "<< elapsed_seconds.count() <<endl;
    }
    if(perm == BEG_PERM){
        e_ans = dd->multiply(e_ans, gg.permuteOperator(n));
    }
    else if(perm == END_PERN){
        e_ans = dd->multiply(gg.permuteOperator(n),e_ans);//this is equivanent to reversing variable order.
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
            gg.lineSet(line, -1, i + j + 1);
            gg.RmatGenerator(Rmat, j+2);
            e_ans = dd->multiply(dd->makeGateDD(Rmat, n, line), e_ans);
            gg.lineReset(line, -1, i + j + 1);
        }
        gg.lineReset(line,i,n - i - 1);//HM: reset line array
    }
    if(perm == END_PERN)
        e_ans = dd->multiply(gg.permuteOperator(n),e_ans);
    
    delete[] line;
    return e_ans;
}
/// QFT in main.
int main(){
    
    if(/* DISABLES CODE */ (false)){//make 'true' to investigate a single QFT (fixed number of bits)
        auto* dd = new dd::Package;
        int n = 3;
        
        //
        //        dd::Edge state3 = dd->makeBasisState(n, 3);
        //        dd::Edge state0 = dd->makeBasisState(n, 0);
        //        dd::Edge state2 = dd->makeBasisState(n, 2);
        //        dd::Edge state1 = dd->makeBasisState(n, 1);
        StateGenerator sg = StateGenerator(dd);
        dd::Edge state = sg.dd_UniformState(n);
        cout << "Input vector with size "<< dd->size(state)<< " \n";
        dd->printVector(state);
        dd->export2Dot(state, "state3before.dot",true, true);
        //  dd->printVector(state);
        //  dd::Edge e_swap =permuteOperator(dd,n);//HM: this can be taken out, here for clarity
        state = dd_QFTV2(dd, n, state,NO_PERM);
        dd->export2Dot(state, "state3after.dot",true, true);
          cout << "Output vector with size "<< dd->size(state)<< " \n";
        dd->printVector(state);
        delete dd;
    }
    
    if(true){//make 'true' to plot the graph (time and node count) as function of varible number.
        int N = 14;//points for graph
        int offset = 2;//starting num of qubits
        float a_et[2][N];//execution time
        int nodecounter[N];
        for (int i = 0; i < N; ++i){
            //Initialize package
            auto* dd = new dd::Package;
            int n = i + offset;//number of bits
            StateGenerator sg = StateGenerator(dd);
            time_point<system_clock> start, end;
            start = system_clock::now();
            dd::Edge e_state = sg.dd_BaseState(n, 2);
//            cout << "Input vector "<<n<<" bits:\n";
//            dd->printVector(e_state);
            e_state = dd_QFTV2(dd, n, e_state, NO_PERM);
//            cout << "Output vector "<<n<<" bits:\n";
//            dd->printVector(e_state);
//            cout << "---\n";
            nodecounter[i] = dd->size(e_state);
            end = system_clock::now();
            duration<float> elapsed_seconds = end - start;
            a_et[0][i] = n;
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
