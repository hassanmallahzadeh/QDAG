//
//  main.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-04.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
//#include "QFT-DDgenerator.hpp"
#include "commonheaders.h"
#include "QFT.hpp"
#include "QFT-Measurement.hpp"
#include "util.h"
#include <stdio.h>
using namespace std::chrono;
using std::cout;
using std::endl;
int main(){
    

    //    short line[2] = {2,-1};
    //    dd::Edge d_test = dd->makeGateDD(Hmat, 2, line);

  //  vector<dd::ComplexValue> vs = {{0.5,0},{0.5,0},{0,0},{0.5,0},{0,0},{0,0},{0.5,0},{0,0}};
//    vector<dd::ComplexValue> vs = {{0,0},{4*0.5,0},{3*0.5,0},{0.5,0},{0.5,0},{0,0},{0,0},{0,0}};
//    //dd::Edge d_test =  gg.dd_Sqrt3State(n);
    if(false){
        auto* dd = new dd::Package;
        int n = 17;
        StateGenerator gg = StateGenerator(dd);
        dd::Edge d_test =  gg.dd_UniformState(n);
        cout<<"input state:"<<endl;
        dd->printVector(d_test,10);
        dd->export2Dot(d_test, "graphs/ d_input2bituniformstate.dot", true, true);
        
        int n_target = {1};
        Measurement mm = Measurement(dd);
        int res = mm.Measure(d_test, n, n_target);
        cout<<"measurement outcome"<<endl;
        cout<< "bit "<< n_target <<": "<< res <<endl;
        cout<<"output state:"<<endl;
        dd->printVector(d_test,10);
        dd->export2Dot(d_test, "graphs/d_output2bituniformstate.dot", true, true);
    }
    
    if(/* DISABLES CODE */ (false)){//make 'true' to investigate a single QFT (fixed number of bits)
        auto* dd = new dd::Package;
        int n = 18;
        QFT qft = QFT(dd);
        StateGenerator sg = StateGenerator(dd);
        dd::Edge state = sg.dd_BaseState(n,0);
        cout<<"input state: \n";
    //    dd->printVector(state);
        time_point<system_clock> start, end;
                  start = system_clock::now();
         
        state =  qft.dd_QFTGNV2(n, state);
        //qft.dd_QFTGN(n, state, NO_PERM);
        end = system_clock::now();
                       duration<float> elapsed_seconds = end - start;
                 
        cout<<"output state: \n";
       // dd->printVector(state);
        cout<< "execution time: "<< elapsed_seconds.count()<<endl;
        //  dd->export2Dot(state, "\DD_Graphs\d_outoutqft2.dot",true, true);
        delete dd;
    }
    
    if(/* DISABLES CODE */ true){//make 'true' to plot the graph (time and node count) as function of varible number.
        int N = 18;//points for graph
        int ntrials = 20;
        int offset = 1;//starting num of qubits
        float nodecounter[N];//type float for averaging
        float a_et[2][N];//dd size, execution time.
        
        for (int j = 0; j < N; ++j){
            nodecounter[j] = 0;
            for (int i = 0; i < 2; ++i){
                a_et[i][j] = 0;
            }
        }
        
        for(int k = 0; k < ntrials; ++k){
            for (int i = 0; i < N; ++i){
                //Initialize package
                auto* dd = new dd::Package;
                int n = i + offset;//number of bits
                StateGenerator sg = StateGenerator(dd);
                dd::Edge state = sg.dd_BaseState(n, pow(2,n) - 1);
                QFT qft = QFT(dd);
                time_point<system_clock> start, end;
                start = system_clock::now();
                state = qft.dd_QFTGNV2(n, state);
                nodecounter[i] += dd->size(state);
                end = system_clock::now();
                duration<float> elapsed_seconds = end - start;
                a_et[0][i] = n;
                a_et[1][i] += elapsed_seconds.count();
                delete dd;
            }
        }
        for(int i = 0; i < N; ++i){
            nodecounter[i] = nodecounter[i] / ntrials;
            a_et[1][i] = a_et[1][i] / ntrials;
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
