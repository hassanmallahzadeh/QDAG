//
//  main.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-04.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
//#include "QFT-DDgenerator.hpp"

#include "QFT.hpp"
#include "util.h"
#include <stdio.h>
using namespace std::chrono;
using std::cout;
using std::endl;
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
        QFT qft = QFT(dd);
        dd::Edge state = sg.dd_UniformState(n);
        cout << "Input vector with size "<< dd->size(state)<< " \n";
        dd->printVector(state);
        dd->export2Dot(state, "state3before.dot",true, true);
        //  dd->printVector(state);
        //  dd::Edge e_swap =permuteOperator(dd,n);//HM: this can be taken out, here for clarity
        state = qft.dd_QFTV2(n, state,NO_PERM);
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
            QFT qft = QFT(dd);
            time_point<system_clock> start, end;
            start = system_clock::now();
            dd::Edge e_state = sg.dd_BaseState(n, 2);
            //            cout << "Input vector "<<n<<" bits:\n";
            //            dd->printVector(e_state);
            e_state = qft.dd_QFTV2(n, e_state, NO_PERM);
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
