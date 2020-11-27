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
#include "shorutil.hpp"
#include "RegisterFactory.hpp"
#include "ProbabilisticPeriodFinder.hpp"
#include "Factorizer.hpp"
#include <chrono>
using namespace std::chrono;
using std::cout;
using std::endl;
void PeriodFinderTest();
void UniformityTest();
void QFTexecutionTimes();
void FinalInRegMeasureWithQFTTest();
int main(){
    lli N = 35;
    Factorizer Alice = Factorizer(N);
    array<lli, 2> f = Alice.Factors();
    printf("Factors of %lli, are %lli and %lli\n", N, f[0], f[1]);
    return 0;
}
void PeriodFinderTest(){
    lli N = 200;
    lli a = 11;
    std::pair<lli,lli> p = {-1,-1};
    int counter = 0;
    do{
        auto* dd = new dd::Package;
        ProbabilisticPeriodFinder pf = ProbabilisticPeriodFinder(N,a,dd);
        p = pf.AttemptReadingMultipleOfInverseOfPeriod();
        if(p.second > 0 && p.first > 0)//0 is not measured.
        {
        ++counter;
        printf("numerator: %lld denominator: %lld run: %d\n",p.first, p.second, counter);
    }
    delete dd;
    } while(!(p.second>0) || shor::modexp(a,p.second,N) != 1/*not good for lli variable to be changed*/);
    printf("period of %lld mod %lld found: %lld in %d run\n",a,N,p.second, counter);
}
void FinalInRegMeasureWithQFTTest(){
    
    lli N = 21;
    lli a = 11;
    int tr = 300;
    vector<int> v;
    int i = 0;
    do{
        auto* dd = new dd::Package;
        ProbabilisticPeriodFinder pf = ProbabilisticPeriodFinder(N,a,dd);
        std::pair<lli,lli> iores = pf.AttemptReadingMultipleOfInverseOfPeriod();
     //   static lli flag = 1;//pick results with same output reg num.
        if(true/*iores.second == flag*/){
            if(v.empty())
                v =  vector<int>(pow(2, pf.ni), 0);//make the info holder vector.
        ++v.at(iores.first);
            ++i;
        }
        delete dd;
    } while(i < tr);
    std::ofstream datafile;
    datafile.open ("/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/QDAG/freqs.txt");
    if(datafile.is_open()){
        for (int i = 0; i < v.size(); ++i){
            datafile << i << " "<< v[i] <<endl;
        }
        datafile.close();
    }
}
/// investigate a single QFT (fixed number of bits)
void UniformityTest(){
    //examine uniformity of probabilities.
    int ntrials = 1024;//
    int n = 4;
    vector<int> vt;//for test of dd_QFTGNV1 when we determine on which qubits QFT is to be performed.
    for(int i = 0; i < n; ++i)
    vt.push_back(i);
    map<string,int> m;
    for(int i = 0; i < ntrials; ++i){
        auto* dd = new dd::Package;
        QFT qft = QFT(dd);
        StateGenerator sg = StateGenerator(dd);
        dd::Edge state = sg.dd_UniformState(n);
        std::random_device rd;
        engine eng(rd());
        qft.dd_QFTGNV1(n, state, NO_PERM, eng,vt);
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
void QFTexecutionTimes(){
    //make 'true' to plot the graph (time and node count) as function of varible number.
    std::random_device rd;
    const int N = 64;//points for graph
    int ntrials = 100;//how many times run simulations for averaging.
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
//            dd::Edge state = sg.dd_BaseState(n,shor::modexp(2, n, std::numeric_limits<lli>::max())-1);
            dd::Edge state = sg.dd_RandomState(n, rd());
         //   dd::Edge state = sg.dd_UniformState(n);
            QFT qft = QFT(dd);
            time_point<system_clock> start, end;
            start = system_clock::now();
            engine eng(rd());
            state = qft.dd_QFTGNV1(n, state, BEG_PERM, eng);
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
