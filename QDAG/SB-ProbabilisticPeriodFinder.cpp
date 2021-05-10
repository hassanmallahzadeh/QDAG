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
SB_PPF::SB_PPF(lli N, lli a, dd::Package* dd): gg(dd) , posControl(1), m_qft(dd){
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
void SB_PPF::phi_Adder(const std::function<int (int)> &acnum, const map<int, bool> &c, short *line, dd::Edge &state){
    for (int i = n; i >= 0; --i){//plus one to accomodate for overflow
        for(int j = i; j >= 0; --j){
            if(acnum(j) == posControl){
                gg.CKRmatGenOrApply(line, m_currentholder(i), i - j + 1, c, nt, &state);
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
void SB_PPF::phi_Subtractor(const std::function<int (int)> &acnum, const map<int, bool> &c, short *line, dd::Edge &state){
    //    Measurement mm = Measurement(dd);
    //     std::random_device rd;
    //     engine eng(rd());
    
    for (int i = 0; i < n + 1; ++i){//plus 1 to accomodate for overflow.
        for(int j = 0; j <= i; ++j){
            if(acnum(j) == posControl){
                gg.CKRmatGenOrApply(line, m_currentholder(i), i - j + 1, c, nt, &state, true);
            }
        }
    }
}
/// Adds in mod N, classical number 'a' to fourier transform of quantum number 'b' in fourier space.
/// Answer will be 'a+b'  mod N in Fourier space. Figure 5 of SB.
/// @param acnum classical number to be added mod N (in base 2) to quantum number b
/// @param currentholder current info holder. b qubits before swap of fig 7 of SB.
/// @param tindex temporary (ancilla) qubit in fig 5 of SB.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. (c1 and c2 in fig 5)
/// @param line 'line' for gate generation
/// @param state input quantum state
void SB_PPF::phi_AdderModN(const std::function<int (int)> &acnum, const std::function<int ()> &tindex, const map<int, bool> &c, short *line, dd::Edge &state){
    gg = GateGenerator (dd);
    vector<int> indice;
    for(int i = 0; i < n + 1; ++i){
        indice.push_back(m_currentholder(i));
    }
    phi_Adder(acnum, c, line, state);
    map<int, bool> temp;//to have an empty map to pass and later use for temp qubit control
    phi_Subtractor(m_base2N, temp, line, state);
    m_qft.dd_QFTV5(nt, state, indice, true);//inverse qft
    gg.CNotGenOrApply(line, tindex(), m_currentholder(n), nt, &state);
    m_qft.dd_QFTV5(nt, state, indice, false);//qft
    //qft
    temp = {{tindex(), true}};
    phi_Adder(m_base2N, temp, line, state);
    phi_Subtractor(acnum, c, line, state);
    m_qft.dd_QFTV5(nt, state, indice, true);//inverse qft
    gg.NotGenOrApply(line, m_currentholder(n), nt, &state);
    gg.CNotGenOrApply(line, tindex(), m_currentholder(n), nt, &state);
    gg.NotGenOrApply(line, m_currentholder(n), nt, &state);
    m_qft.dd_QFTV5(nt, state, indice, false);//qft
    phi_Adder(acnum, c, line, state);
}
/// Subtracts in mod N, classical number 'a' from fourier transform of quantum number 'b' in fourier space.
/// Answer will be 'b - a'  mod N in Fourier space. Figure 5 of SB (for addition case).
/// @param acnum classical number to be subtracted mod N (in base 2) from quantum number b
/// @param bindice b indice generated by a lambda. n + 1 th element is the overflow qubit.
/// @param tindex temporary (ancilla) qubit in fig 5 of SB.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. (c1 and c2 in fig 5)
/// @param line 'line' for gate generation
/// @param state input quantum state
void SB_PPF::phi_SubtractorModN(const std::function<int (int)> &acnum, const std::function<int ()> &tindex, const map<int, bool> &c, short *line, dd::Edge &state){
    gg = GateGenerator (dd);
    vector<int> indice;
    for(int i = 0; i < n + 1; ++i){
        indice.push_back(m_currentholder(i));
    }
    phi_Subtractor(acnum, c, line, state);
    m_qft.dd_QFTV5(nt, state, indice, true);//qftreverse
    gg.NotGenOrApply(line, m_currentholder(n), nt, &state);
    gg.CNotGenOrApply(line, tindex(), m_currentholder(n), nt, &state);
    gg.NotGenOrApply(line, m_currentholder(n), nt, &state);
    m_qft.dd_QFTV5(nt, state, indice, false);//inverse qft
    phi_Adder(acnum, c, line, state);
    map<int, bool> temp;//temp control
    temp = {{tindex(), true}};
    phi_Subtractor(m_base2N, temp, line, state);
    m_qft.dd_QFTV5(nt, state, indice, true);
    gg.CNotGenOrApply(line, tindex(), m_currentholder(n), nt, &state);
    m_qft.dd_QFTV5(nt, state, indice, false);//inverse qft
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
void SB_PPF::phi_CMultiplier(const lli &a0cnum, const std::function<int ()> &tindex, map<int,bool> c, short *line, dd::Edge &state){
    vector<int> indice;
    int shiftedn = n + 1;
    for(int i = 0; i < shiftedn; ++i){
        indice.push_back(m_currentholder(i));
    }
    static bool test = false;
    if(!test){
        dd->export2Dot(state, "before.dot");
    }
    m_qft.dd_QFTV5(nt, state, indice, false);//qft
    if(!test){
        dd->export2Dot(state, "after.dot");
        test = true;
    }
    for(int i = 0; i < n; ++i){//fig6 of paper, repeated shift and apply
        //TODO: make next four opertations not do redundant calculations.
        lli tempa0cnum = a0cnum;
        tempa0cnum = tempa0cnum << i;
        tempa0cnum = tempa0cnum % N;
        vector<bool> a0clb2 = shor::base2rep(tempa0cnum, shiftedn);//classical value in base 2, n digits.
        auto a0cnumbase2 = [&a0clb2](int i){return ((a0clb2.empty() || i >= a0clb2.size()) ? 0 : a0clb2[i]); };//return cnum's digits in base 2.
        c.insert({m_currentholder(i),posControl});//watch for map random insert.
        phi_AdderModN(a0cnumbase2, tindex ,c , line, state);
        c.erase(--c.end());
    }
    m_qft.dd_QFTV5(nt, state, indice, true);//inverse qft
}
/// Divides quantum number x by classical number 'a' and subtracts the result from quantum number in 'b' register, if 'c'  qubit is 1, else leave b qubits unchanged.
/// Figure 6 of SB. x is unchanged in any case
/// @param a0cnum classical number to divide x, mod N (in base 2)
/// @param bindice b indice generated by a lambda.  n + 1 th element is the overflow qubit.
/// @param xindice x indice generated by a lambda.  this is the qunatum number to be divided by classical number a0cnum
/// @param tindex temporary (ancilla) qubit in fig 5 of SB.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. map leaves potential for more than 1 control qubit.
/// @param line 'line' for gate generation
/// @param state input quantum state
void SB_PPF::phi_CDivider(const lli &acnum, const std::function<int ()> &tindex, map<int,bool> c, short *line, dd::Edge &state){
    vector<int> indice;
    for(int i = 0; i < n + 1; ++i){
        indice.push_back(m_currentholder(i));
    }
    m_qft.dd_QFTV5(nt, state, indice, true);//inverse qft
    for(int i = n - 1 ; i >= 0; --i){//fig5 of paper, repeated shift and apply
        //TODO: make next four opertations not do redundant calculations.
        lli tempa0cnum = acnum;
        tempa0cnum = tempa0cnum << i;
        tempa0cnum = tempa0cnum % N;
        vector<bool> a0clb2 = shor::base2rep(tempa0cnum, n);//classical value in base 2, n digits.
        auto a0cnumbase2 = [&a0clb2](int i){return ((a0clb2.empty() || i >= a0clb2.size()) ? 0 : a0clb2[i]); };//return cnum's digits in base 2
        c.insert({m_currentholder(i),posControl});//watch for map random insert.
        phi_SubtractorModN(a0cnumbase2,tindex ,c , line, state);//temp remove
        c.erase(--c.end());
    }
    m_qft.dd_QFTV5(nt, state, indice, false);//qft
}
/// 'Multiplies quantum number x to classical number 'a' and adds the result to quantum number in 'b' register (|0> in this case), puts result to x register, resets b to inital value' if 'c'  qubit is 1, else leave x and b qubits unchanged.
/// Figure 7 of SB and equation 3.
/// @param acnum classical number to be multiplied mod N (in base 2)
/// @param currentholder current info holder. b qubits before swap of fig 7 of SB.
/// @param tindex temporary (ancilla) qubit in fig 5 of SB.
/// @param c map of control qubit indices. value 'true' measns positive control and value 'false' means negative control. map leaves potential for more than 1 control qubit.
/// 'c' qubits in Fig 7 of SB
/// @param line 'line' for gate generation
/// @param state input quantum state
   void SB_PPF::phi_CUa(const lli &acnum, const std::function<int ()> &tindex, map<int,bool> c, short *line, dd::Edge &state){
  
    phi_CMultiplier(acnum,tindex, c, line, state);

    lli invacnum = shor::modInverse(acnum, N);
       m_bbased = !m_bbased;//swap
    phi_CDivider(invacnum, tindex, c, line, state);//temp, change
    
    //  phi_CMultiplier(acnum,bindice,xindice,tindex, c, line, state);
}

std::pair<lli,lli> SB_PPF::AttemptReadingMultipleOfInverseOfPeriod(){
    int m =  (int)shor::base2rep(N*N).size();
    short* line = new short[nt];
    gg = GateGenerator (dd);
    gg.lineClear(line, nt);
    Measurement mm = Measurement(dd);
    dd::Edge state = StateGenerator(dd).dd_BaseState(nt, 0);//start by setting all qubits to zero.
    std::function<int (void)> cindex =  [](){return 0;};//top control qubit on which measurements are performed to read and report the result. Fig 8 of SB.
    std::function<int (void)> tindex =  [](){return 1;};//temporary (ancilla) qubit in fig 5 of SB.
    auto b0indice = [n = this->n](int i){assert(i >= 0 && i < n);
        return i + 2;}; //qft is applied to this registers.
    auto overflow = [n = this->n](){return n + 2;};//by default over flow set one above b indices
    auto xindice = [n = this->n](int i){return n + i + 3;}; // x qubits used in multiplier and U_a. They store result (a*x Mod N) in U_a. n qubits.
    m_currentholder = [n = this->n, m_bbased= &this->m_bbased, overflow, xindice, b0indice](int i){//x and b qubits are swapped as per fig 7 of SB. But, we must fix one qubit for overflow on current qubit number holder.
        if(i < n){
            if(m_bbased){
                return b0indice(i);
            }
            else{
                return xindice(i);
            }
        }
        else if(i == n){
            return overflow();
        }
        else {
            assert(0);
        }
    };
    m_bbased = false;
    gg.NotGenOrApply(line, m_currentholder(0), nt, &state);//put 1 in x register to start(fig8 of SB)
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
        if(i % 2 == 0){
            assert(m_bbased == true);
        }
        else{//swap
            assert(m_bbased == false);
        }
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
        //  dd->export2Dot(state, "test1.dot");
        
        int mres = mm.Measure(state, nt, cindex(), eng);//warning: watch for random number generator something fishy might be here. I once saw something...
        //   dd->export2Dot(state, "test2.dot");
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
int main(){
    auto* dd = new dd::Package;
    SB_PPF pf = SB_PPF(3,2,dd);
    pf.AttemptReadingMultipleOfInverseOfPeriod();
    delete dd;
}
