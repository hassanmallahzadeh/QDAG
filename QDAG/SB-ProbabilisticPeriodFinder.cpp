//
//  SB-ProbabilisticPeriodFinder.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2021-03-18.
//  Copyright © 2021 Hassan Mallahzadeh. All rights reserved.
// This file implements Shor's algorithm based on Stephane Beauregard's paper "Circuit for Shor's algorithm using 2n+3 qubits". Paper refered to 'SB' below.
#include "SB-ProbabilisticPeriodFinder.hpp"
#include "QFT-DDgenerator.hpp"
#include "shorutil.hpp"

/// Constructor for SB Shor's algorithm scheme.
/// @param N Number to be find period of base a
/// @param a base for finding the period on N
/// @param dd pointer to qmdd object.
SB_ProbabilisticPeriodFinder::SB_ProbabilisticPeriodFinder(lli N, lli a, dd::Package* dd): gg(dd) , posControl(1), m_qft(dd){
    dd ? this->dd = dd : dd = new dd::Package;
    this->N = N;
    this->a = a;
    auto base2N = shor::base2rep(N);//base 2 representation of N
    n = static_cast<int>(base2N.size());//num qubits needed to represent N.
    nt = 2 * n + 3;//n for x register. +1 temporary for modNadder. + (n + 1) for auxiliary '0' qubits in multiplier. +1 for control qubit of U-a gates.
    
    m_base2N = [base2N](int i){return base2N.empty() ? 0 : base2N[i]; };
    //return cnum's digits in base 2.
    //wrapper to read cnum classic bits.
}
/// Adds classical number 'a' to fourier transform of quantum number 'b' in fourier space. Fig 2 & 3 of SB.
/// Answer will be 'a+b' in Fourier space. Figure 3 (3 implemented, a is classic number) of SB.
/// @param acnum a0, classic number value in base 2.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. This is needed for modN adder (fig 5)
/// @param line 'line' for gate generation
/// @param state input quantum state
void SB_ProbabilisticPeriodFinder::phi_Adder(const std::function<int (int)> &acnum, const map<int, bool> &c, short *line, dd::Edge &state){
    for (int i = n; i >= 0; --i){//plus one to accomodate for overflow
        for(int j = i; j >= 0; --j){
            if(acnum(j) == posControl){
                gg.CKRmatGenOrApply(line, currentHolder(i), i - j + 1, c, nt, &state);
            }
        }
    }
}

/// Subtracts classical number 'a' from fourier transform of quantum number 'b' in fourier space. Fig 2 & 3 of SB, run in reverse.
/// Answer will be 'b - a' in Fourier space. Figure 3 (3 implemented, a is classic number) of SB, run in reverse.
//// @param acnum a0, classic number value in base 2.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. This is needed for modN adder (fig 5)
/// @param line 'line' for gate generation
/// @param state input quantum state
void SB_ProbabilisticPeriodFinder::phi_Subtractor(const std::function<int (int)> &acnum, const map<int, bool> &c, short *line, dd::Edge &state){
    //    Measurement mm = Measurement(dd);
    //     std::random_device rd;
    //     engine eng(rd());
    
    for (int i = 0; i < n + 1; ++i){//plus 1 to accomodate for overflow.
        for(int j = 0; j <= i; ++j){
            if(acnum(j) == posControl){
                gg.CKRmatGenOrApply(line, currentHolder(i), i - j + 1, c, nt, &state, true);
            }
        }
    }
}
/// Adds in mod N, classical number 'a' to fourier transform of quantum number 'b' in fourier space.
/// Answer will be 'a+b'  mod N in Fourier space. Figure 5 of SB.
/// @param acnum classical number to be added mod N (in base 2) to quantum number b
/// @param tindex temporary (ancilla) qubit in fig 5 of SB.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. (c1 and c2 in fig 5)
/// @param line 'line' for gate generation
/// @param state input quantum state
void SB_ProbabilisticPeriodFinder::phi_AdderModN(const std::function<int (int)> &acnum, const std::function<int ()> &tindex, const map<int, bool> &c, short *line, dd::Edge &state){
    gg = GateGenerator (dd);
    phi_Adder(acnum, c, line, state);
    map<int, bool> temp;//to have an empty map to pass and later use for temp qubit control
    phi_Subtractor(m_base2N, temp, line, state);
    m_qft.dd_QFTV5(nt, state, m_qftindices, true);//inverse qft
    gg.CNotGenOrApply(line, tindex(), currentHolder(n), nt, &state);
    m_qft.dd_QFTV5(nt, state, m_qftindices, false);//qft
    //qft
    temp = {{tindex(), true}};
    phi_Adder(m_base2N, temp, line, state);
    phi_Subtractor(acnum, c, line, state);
    m_qft.dd_QFTV5(nt, state, m_qftindices, true);//inverse qft
    gg.NotGenOrApply(line, currentHolder(n), nt, &state);
    gg.CNotGenOrApply(line, tindex(), currentHolder(n), nt, &state);
    gg.NotGenOrApply(line, currentHolder(n), nt, &state);
    m_qft.dd_QFTV5(nt, state, m_qftindices, false);//qft
    phi_Adder(acnum, c, line, state);
}
/// Subtracts in mod N, classical number 'a' from fourier transform of quantum number 'b' in fourier space.
/// Answer will be 'b - a'  mod N in Fourier space. Figure 5 of SB (for addition case).
/// @param acnum classical number to be subtracted mod N (in base 2) from quantum number b
/// @param tindex temporary (ancilla) qubit in fig 5 of SB.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. (c1 and c2 in fig 5)
/// @param line 'line' for gate generation
/// @param state input quantum state
void SB_ProbabilisticPeriodFinder::phi_SubtractorModN(const std::function<int (int)> &acnum, const std::function<int ()> &tindex, const map<int, bool> &c, short *line, dd::Edge &state){
    gg = GateGenerator (dd);
    phi_Subtractor(acnum, c, line, state);
    m_qft.dd_QFTV5(nt, state, m_qftindices, true);//qftreverse
    gg.NotGenOrApply(line, currentHolder(n), nt, &state);
    gg.CNotGenOrApply(line, tindex(), currentHolder(n), nt, &state);
    gg.NotGenOrApply(line, currentHolder(n), nt, &state);
    m_qft.dd_QFTV5(nt, state, m_qftindices, false);//inverse qft
    map<int, bool> temp;//temp control
    temp = {{tindex(), true}};
    phi_Adder(acnum, c, line, state);
    phi_Subtractor(m_base2N, temp, line, state);
    m_qft.dd_QFTV5(nt, state, m_qftindices, true);
    gg.CNotGenOrApply(line, tindex(), currentHolder(n), nt, &state);
    m_qft.dd_QFTV5(nt, state, m_qftindices, false);//inverse qft
    temp.clear();
    phi_Adder(m_base2N, temp, line, state);
    phi_Subtractor(acnum, c, line, state);//
}
/// Multiplies quantum number x to classical number 'a' and adds the result to quantum number in 'b' register,if 'c'  qubit is 1, else leave b qubits unchanged.
/// Figure 6 of SB.. x is unchanged in any case
/// @param a0cnum classical number to be multiplied by x mod N (in base 2)
/// @param bindice b indice generated by a lambda.  n + 1 th element is the overflow qubit.
/// @param xindice x indice generated by a lambda.  this is the qunatum number to be multiplied by classical number a0cnum
/// @param tindex temporary (ancilla) qubit in fig 5 of SB.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. map leaves potential for more than 1 control qubit.
/// @param line 'line' for gate generation
/// @param state input quantum state
void SB_ProbabilisticPeriodFinder::phi_CMultiplier(const lli &a0cnum, const std::function<int ()> &tindex, map<int,bool> c, short *line, dd::Edge &state){
    
    int shiftedn = n + 1;
    m_qftindices.clear();
    for(int i = 0; i < shiftedn; ++i){
        m_qftindices.push_back(currentHolder(i));
    }
    m_qft.dd_QFTV5(nt, state, m_qftindices, false);//qft
    for(int i = 0; i < n; ++i){//fig6 of paper, repeated shift and apply
        //TODO: make next four opertations not do redundant calculations.
        lli tempa0cnum = a0cnum;
        tempa0cnum = tempa0cnum << i;
        tempa0cnum = tempa0cnum % N;
        vector<bool> a0clb2 = shor::base2rep(tempa0cnum, shiftedn);//classical value in base 2, n digits.
        auto a0cnumbase2 = [&a0clb2](int i){return ((a0clb2.empty() || i >= a0clb2.size()) ? 0 : a0clb2[i]); };//return cnum's digits in base 2.
        m_bbased = !m_bbased;//when b is current data holder, x gives control and vice versa.
        int currentholder = currentHolder(i);
        c.insert({currentholder,posControl});
        m_bbased = !m_bbased;
        phi_AdderModN(a0cnumbase2, tindex ,c , line, state);
        c.erase(c.find(currentholder));
    }
    m_qft.dd_QFTV5(nt, state, m_qftindices, true);//inverse qft
    gg.NotGenOrApply(line, c.begin()->first, nt, &state);
    for(int i = 0; i < n; ++i){
        m_bbased = !m_bbased;//when b is current data holder, x gives control and vice versa.
        int currentcontrol = currentHolder(i);
        m_bbased = !m_bbased;//when b is current data holder, x gives control and vice versa.
        gg.ToffoliGenOrApply(line, currentHolder(i), c.begin()->first, currentcontrol, nt, &state);
    }
    gg.NotGenOrApply(line, c.begin()->first, nt, &state);
}
/// Divides quantum number x by classical number 'a' and subtracts the result from quantum number in 'b' register, if 'c'  qubit is 1, else leave b qubits unchanged.
/// Figure 6 of SB. x is unchanged in any case
/// @param tindex temporary (ancilla) qubit in fig 5 of SB.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. map leaves potential for more than 1 control qubit.
/// @param line 'line' for gate generation
/// @param state input quantum state
void SB_ProbabilisticPeriodFinder::phi_CDivider(const lli &acnum, const std::function<int ()> &tindex, map<int,bool> c, short *line, dd::Edge &state){
    gg.NotGenOrApply(line, c.begin()->first, nt, &state);
    for(int i = n - 1; i >= 0; --i){
        m_bbased = !m_bbased;//when b is current data holder, x gives control and vice versa.
        int currentcontrol = currentHolder(i);
        m_bbased = !m_bbased;//when b is current data holder, x gives control and vice versa.
        gg.ToffoliGenOrApply(line, currentHolder(i), c.begin()->first, currentcontrol, nt, &state);
    }
    gg.NotGenOrApply(line, c.begin()->first, nt, &state);
    int shiftedn = n + 1;

    m_qftindices.clear();
    for(int i = 0; i < shiftedn; ++i){
        m_qftindices.push_back(currentHolder(i));
    }
    m_qft.dd_QFTV5(nt, state, m_qftindices, false);//
    for(int i = n - 1 ; i >= 0; --i){//fig5 of paper, repeated shift and apply
        //TODO: make next four opertations not do redundant calculations.
        lli tempa0cnum = acnum;
        tempa0cnum = tempa0cnum << i;
        tempa0cnum = tempa0cnum % N;
        vector<bool> a0clb2 = shor::base2rep(tempa0cnum, shiftedn);//classical value in base 2, n digits.
        auto a0cnumbase2 = [&a0clb2](int i){return ((a0clb2.empty() || i >= a0clb2.size()) ? 0 : a0clb2[i]); };//return cnum's digits in base 2
        m_bbased = !m_bbased;//when b is current data holder, x gives control and vice versa.
        int currentholder = currentHolder(i);
        c.insert({currentholder,posControl});
         m_bbased = !m_bbased;
        phi_SubtractorModN(a0cnumbase2,tindex ,c , line, state);
        c.erase(c.find(currentholder));
    }
    m_qft.dd_QFTV5(nt, state, m_qftindices, true);//inverse qft
}
/// 'Multiplies quantum number x to classical number 'a' and adds the result to quantum number in 'b' register (|0> in this case), puts result to x register, resets b to inital value' if 'c'  qubit is 1, else leave x and b qubits unchanged.
/// Figure 7 of SB and equation 3.
/// @param acnum classical number to be multiplied mod N (in base 2)
/// @param tindex temporary (ancilla) qubit in fig 5 of SB.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. map leaves potential for more than 1 control qubit.
/// 'c' qubits in Fig 7 of SB
/// @param line 'line' for gate generation
/// @param state input quantum state
   void SB_ProbabilisticPeriodFinder::phi_CUa(const lli &acnum, const std::function<int ()> &tindex, map<int,bool> c, short *line, dd::Edge &state){
    phi_CMultiplier(acnum,tindex, c, line, state);
    lli invacnum = shor::modInverse(acnum, N);
       m_bbased = !m_bbased;
//       vector<int> vec0;
//       vector<int> vec1;
//       for(int i = 0 ; i < n ; ++i)
//      vec0.push_back(currentHolder(i));//watch for map random insert.
//      m_bbased = !m_bbased;//when b is current data holder, x gives control and vice versa.
//       for(int i = 0 ; i < n ; ++i)
//      vec1.push_back(currentHolder(i));//watch for map random insert.
//       m_bbased = !m_bbased;
//       state = gg.swapRegistersOnState(nt,vec0,vec1,state,c);
    phi_CDivider(invacnum, tindex, c, line, state);//temp, change
    
    //  phi_CMultiplier(acnum,bindice,xindice,tindex, c, line, state);
}

std::pair<lli,lli> SB_ProbabilisticPeriodFinder::AttemptReadingMultipleOfInverseOfPeriod(){
    int m =  (int)shor::base2rep(N*N).size();
    short* line = new short[nt];
    gg = GateGenerator (dd);
    gg.lineClear(line, nt);
    Measurement mm = Measurement(dd);
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    std::function<int (void)> cindex =  [=](){return nt - 1;};//top control qubit on which measurements are performed to read and report the result. Fig 8 of SB.
    std::function<int (void)> tindex =  [](){return 0;};//temporary (ancilla) qubit in fig 5 of SB.
    m_bbased = false;
    gg.NotGenOrApply(line, currentHolder(0), nt, &state);//put 1 in x register to start(fig8 of SB)
    m_bbased = true;
    lli cnum = a;//classical number to be multiplied by quantum number.
    std::random_device rd;
    engine eng(rd());
    vector<int> vmres;//binary representation of one full set of measurements
    vector<lli> factors;//stores classical multiplication factor (a^2^i) , fig 8 of SB
    map<int,bool> c = {{cindex(), true}};//this map mechanism is kind of silly but lets not worry about that.
    for (int i = 0; i < m; ++i){
        factors.push_back(cnum);
        cnum = cnum * cnum % N;
    }
    for(int i = 0; i < m; ++i){
        int ii = m - i - 1;
        //fig 1 of HRS.
        gg.HadGenOrApply(line, cindex(), nt, &state);
        phi_CUa(factors[ii], tindex, c, line, state);
        dd::Matrix2x2 Rmat;// put rotation gates in place
        gg.lineSet(line, cindex(), -1);
        for(int j = 0; j < i; ++j){
            if(vmres[j]){
                gg.RmatGenerator(Rmat, i - j + 1);
                state = dd->multiply(dd->makeGateDD(Rmat, nt, line), state);
            }
        }
        gg.lineReset(line, cindex(), -1);
        gg.HadGenOrApply(line, cindex(), nt, &state);
        int mres = mm.Measure(state, nt, cindex(), eng);//warning: watch for random number generator something fishy might be here. I once saw something...
        vmres.push_back(mres);
        if(mres == posControl){
            gg.NotGenOrApply(line, cindex(), nt, &state);
        }
    }
    lli res = shor::base2to10(vmres, false);
    std::pair<lli,lli> p = shor::contfrac(res, m,n);
    return p;
}
///// 'out_Circuit, top circuit of Shor's algorithm on SB scheme, Fig 8.
///// 'Multiplies quantum number x to classical number 'a' and adds the result to quantum number in 'b' register (|0> in this case), puts result to x register, resets b to inital value' if 'c'  qubit is 1, else leave x and b qubits unchanged.
///// Figure 7 of SB and equation 3.
///// @param bindice b indice generated by a lambda (|0> in this case) result of multiplication of x and a are put in b registers.
///// @param tindex temporary (ancilla) qubit in fig 5 of SB.
///// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. map leaves potential for more than 1 control qubit.
///// 'c' qubits in Fig 7 of SB
///// @param line 'line' for gate generation
///// @param state input quantum state
///// @param xindice  quantum register to be multiplied by acnum mod N. x indice (for input quantum number)
//void SB_PPF::out_Circuit(const std::function<int (int)> &bindice, const std::function<int (int)> &xindice, const std::function<int ()> &tindex, const std::function<int ()> &mindex, map<int,bool> c, short *line, dd::Edge &state){
//    gg.NotGenOrApply(line, xindice(0), nt, &state);//put 1 in x register (fig 8 of SB)
//    lli cnum = a;
//    vector<int> vmres;//binary representation of one full set of measurements
//    vector<lli> factors;//stores classical multiplication factor, fig 6 of VBE
//    int m =  (int)shor::base2rep(N*N).size();
//      for (int i = 0; i < m; ++i){
//        factors.push_back(cnum);
//        cnum = cnum * cnum % N;
//        gg.HadGenOrApply(line, mindex(), nt, &state);
//          int ii = m - i - 1;
//          if(ii % 2 == 0){
//          
//      }
//          else{
//               }
//    }
//    
//    
//   //fig 1 of HRS.
//
//    gg.HadGenOrApply(line, xindex(), nt, &state);
//    if(ii % 2 == 0){
//    }
//   CMultiplierModuloNHalfClassic(a0Nbase2, b0indice, c0indice, c1indice, factors[ii], line, xindex(), nt, state, tindex, xpindice);
//   lli invcnum = shor::modInverse(factors[ii], N);
//   InvCMultiplierModuloNHalfClassic(a0Nbase2, xpindice, c0indice, c1indice, invcnum, line, xindex(), nt, state, tindex, b0indice);
//}
//    
//    for(int i = 0 ; i < 2 * n; ++i){//watch for number of qubits in this loop. we are using one extra for overflow then doubling.
//        cnum *= cnum;
//        phi_CUa(cnum, bindice, xindice, tindex, c, line, state);
//    }
//    phi_CMultiplier(cnum,bindice,xindice,tindex,c, line, state);
//    lli invacnum = shor::modInverse(cnum, N);
//    phi_CDivider(invacnum, xindice, bindice, tindex, c, line, state);
//}
int SB_ProbabilisticPeriodFinder::currentHolder(int i){
    if(i < n){
        if(
           m_bbased){
            return i + 1;//b qubits
        }
        else{
            return n + i + 2;//x qubits
        }
    }
    else if(i == n){
        return n + 1;//overflow qubit
    }
    else {
        assert(0);
        return -1;
    }
}
//int main(){
//    auto* dd = new dd::Package;
//    SB_PPF pf = SB_PPF(3,2,dd);
//    pf.AttemptReadingMultipleOfInverseOfPeriod();
//    delete dd;
//}
