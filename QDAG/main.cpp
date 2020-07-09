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
          
    if(/* DISABLES CODE */ (true)){//make 'true' to investigate a single QFT (fixed number of bits)
        auto* dd = new dd::Package;
        int n = 3;        
        QFT qft = QFT(dd);
        dd::Edge d_qft =  qft.dd_QFTV3(n, END_PERN);
        dd->export2Dot(d_qft, "qftinverse3.dot",false, true);
        delete dd;
    }
    
    if(/* DISABLES CODE */ false){//make 'true' to plot the graph (time and node count) as function of varible number.
        int N = 7;//points for graph
        int offset = 2;//starting num of qubits
        float a_et[2][N];//execution time
        int nodecounter[N];
        for (int i = 0; i < N; ++i){
            //Initialize package
            auto* dd = new dd::Package;            
            int n = i + offset;//number of bits
            QFT qft = QFT(dd);
            time_point<system_clock> start, end;
            start = system_clock::now();
            dd::Edge e_qft = qft.dd_QFTV1(n, BEG_PERM);
            nodecounter[i] = dd->size(e_qft);
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
