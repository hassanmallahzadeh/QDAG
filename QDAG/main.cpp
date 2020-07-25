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
#include <random>
using namespace std::chrono;
using std::cout;
using std::endl;
int main(){
    if(false){
        auto* dd = new dd::Package;
        int n = 2;
        StateGenerator sg = StateGenerator(dd);
        dd::Edge state = sg.dd_BaseState(n,2);
        cout<<"input state: \n";
        dd->printVector(state);

        GateGenerator  gg = GateGenerator(dd);
        dd::Edge perm_gate = gg.permuteOperator(n);
        state = dd->multiply(perm_gate, state);
        cout<<"output state: \n";
        dd->printVector(state);

        dd->export2Dot(perm_gate, "d_permgate4.dot", false, true);
    }

    if(/* DISABLES CODE */ (true)){//make 'true' to investigate a single QFT (fixed number of bits)
        int ntrials = 8 * 1024;//
        int n = 3;
        map<string,int> m;
        for(int i = 0; i < ntrials; ++i){
            auto* dd = new dd::Package;
            QFT qft = QFT(dd);
            StateGenerator sg = StateGenerator(dd);
            dd::Edge state = sg.dd_BaseState(n, 0);
            std::random_device rd;
              std::mt19937 engine(rd());
            state =  qft.dd_QFTGNV1(n, state, BEG_PERM, engine);
            dd::NodePtr p = state.p;
            string s = "";
            while(p != dd->terminalNode){
                for(int i = 0; i < dd::NEDGE; ++i){
                    if(i%dd::RADIX!=0)
                        continue;
                    if(!dd->cn.equalsZero(p->e[i].w)){
                        p = p->e[i].p;
                        s = s + std::to_string(i/dd::RADIX);
                        break;
                    }
                }
            }

            if(m.find(s) == m.cend()){
                m.insert(make_pair(s,1));
            }
            else{
                m[s] = m[s] + 1;
            }
        delete dd;
        }
          for(auto it = m.cbegin(); it != m.cend(); ++it){
            cout<<it->first <<": "<<it->second<<endl;
        }
    }

    if(/* DISABLES CODE */ true){//make 'true' to plot the graph (time and node count) as function of varible number.
        int N = 30;//points for graph
        int ntrials = 1;//how many times run simulations for averaging.
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
                std::random_device rd;
                dd::Edge state = sg.dd_UniformState(n);
                QFT qft = QFT(dd);
                time_point<system_clock> start, end;
                start = system_clock::now();
                state = qft.dd_QFTV2(n, state, NO_PERM);
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

